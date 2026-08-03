#pragma once

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/mechanical/mechanical_parameters_2d.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
class _temperature_condition final {
    const mesh::mesh_2d<T>& _mesh;

    explicit _temperature_condition(const mesh::mesh_2d<T>& mesh) : _mesh{mesh} {}

    template<class Hooke, class Thermal_Strain>
    std::array<T, 2> operator()(const Hooke& hooke_matrix, 
                                const Thermal_Strain& thermal_strain,
                                const size_t e, const size_t i) const {
        using namespace metamath::operators;
        std::array<T, 2> integral = {};
        size_t qshift = _mesh.quad_shift(e);
        const auto& el = _mesh.container().element_2d(e);
        for(const size_t q : el.qnodes()) {
            const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshift] : 
                                                       std::get<Constant>(hooke_matrix);
            const auto thermal_stress = calc_stress<T>(hooke, thermal_strain[qshift]);
            const auto wdN = el.weight(q) * _mesh.derivatives(e, i, q);
            integral[X] += wdN[X] * thermal_stress[XX] + wdN[Y] * thermal_stress[XY];
            integral[Y] += wdN[Y] * thermal_stress[YY] + wdN[X] * thermal_stress[XY];
            ++qshift;
        }
        return integral;
    }

    template<class Hooke, class Thermal_Strain>
    std::array<T, 2> operator()(const Hooke& hooke_matrix, const Thermal_Strain& thermal_strain,
                                const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
                                const size_t eL, const size_t eNL, const size_t iL) const {
        using namespace metamath::operators;
        std::array<T, 2> integral = {};
        const auto& elL = _mesh.container().element_2d(eL);
        const auto& elNL = _mesh.container().element_2d(eNL);
        size_t qshiftL = _mesh.quad_shift(eL);
        for(const size_t qL : elL.qnodes()) {
            std::array<T, 3> inner_integral = {};
            size_t qshiftNL = _mesh.quad_shift(eNL);
            const std::array<T, 2>& qcoordL = _mesh.quad_coord(eL, qL);
            for(const size_t qNL : elNL.qnodes()) {
                const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshiftNL] : 
                                                           std::get<Constant>(hooke_matrix);
                const T weight = elNL.weight(qNL) * influence(qcoordL, _mesh.quad_coord(qshiftNL)) * _mesh.jacobian(qshiftNL);
                inner_integral += weight * calc_stress<T>(hooke, thermal_strain[qshiftNL]);
                ++qshiftNL;
            }
            const auto wdN = elL.weight(qL) * _mesh.derivatives(eL, iL, qL);
            integral[X] += wdN[X] * inner_integral[XX] + wdN[Y] * inner_integral[XY];
            integral[Y] += wdN[Y] * inner_integral[YY] + wdN[X] * inner_integral[XY];
            ++qshiftL;
        }
        return integral;
    }

public:
    template<std::floating_point U>
    friend void temperature_condition(std::vector<std::array<U, 2>>& f, const mesh::mesh_2d<U>& mesh,
                                      const evaluated_mechanical_parameters<U>& parameters);
};

template<std::floating_point T>
void temperature_condition(std::vector<std::array<T, 2>>& f, const mesh::mesh_2d<T>& mesh,
                           const evaluated_mechanical_parameters<T>& parameters) {
    const _temperature_condition<T> integrator{mesh};
    const auto process_node = mesh.process_nodes();
#pragma omp parallel for default(none) shared(f, mesh, parameters, integrator, process_node) schedule(dynamic)
    for(size_t node = process_node.front(); node < *process_node.end(); ++node) {
        std::array<T, 2> integral = {};
        for(const size_t eL : mesh.elements(node)) {
            const auto& group = mesh.container().group(eL);
            const auto& [model, physical] = parameters.at(group);
            std::visit(metamath::types::visitor{
                [](const auto&, const std::monostate) {},
                [&](const auto& hooke, const auto& thermal_strain) {
                    using namespace metamath::operators;
                    const size_t iL = mesh.global_to_local(eL, node);
                    if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                        for(const size_t eNL : mesh.neighbours(eL))
                            integral += integrator(hooke, thermal_strain, model.influence, eL, eNL, iL);
                        integral *= nonlocal::nonlocal_weight(model.local_weight);
                    }
                    integral += model.local_weight * integrator(hooke, thermal_strain, eL, iL);
                }
            }, physical.elastic, physical.thermal_strain);
        }
        using namespace metamath::operators;
        f[node] += integral;
    }
}

}
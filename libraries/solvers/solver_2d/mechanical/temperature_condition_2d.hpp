    #pragma once

#include <mesh/mesh_2d/mesh_2d_utils.hpp>

#include <Eigen/Dense>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I>
class _temperature_condition final {
    const std::vector<T> _delta_temperature;
    const mesh::mesh_2d<T, I>& _mesh;

    explicit _temperature_condition(const mesh::mesh_2d<T, I>& mesh,
                                    const std::vector<T>& delta_temperature)
        : _delta_temperature{nonlocal::mesh::utils::nodes_to_qnodes(mesh, delta_temperature)}
        , _mesh{mesh} {}

    template<class Hooke, class Expansion>
    std::array<T, 3> calc_stress(const Hooke& hooke_matrix, const Expansion& thermal_expansion, const size_t qshift) const {
        using namespace metamath::functions;
        using IHM = evaluated_isotropic_hook_matrix_t<T>;
        using OHM = evaluated_orthotropic_hook_matrix_t<T>;
        using AHM = evaluated_anisotropic_hook_matrix_t<T>;
        using ITE = evaluated_isotropic_thermal_expansion_t<T>;
        using OTE = evaluated_orthotropic_thermal_expansion_t<T>;
        using ATE = evaluated_anisotropic_thermal_expansion_t<T>;
        const auto& hooke = hooke_matrix.index() ? std::get<Variable>(hooke_matrix)[qshift] : 
                                                   std::get<Constant>(hooke_matrix);
        const auto& expansion = thermal_expansion.index() ? std::get<Variable>(thermal_expansion)[qshift] : 
                                                            std::get<Constant>(thermal_expansion);
        const auto thermal_strain = _delta_temperature[qshift] * expansion;
        if constexpr (std::is_same_v<Hooke, IHM>) {
            using namespace isotropic_indices;
            if constexpr (std::is_same_v<Expansion, ITE>) {
                const auto thermal_stress = (hooke[_11] + hooke[_12]) * thermal_strain;
                return {thermal_stress, thermal_stress};
            } else if constexpr (std::is_same_v<Expansion, OTE>) {
                return {hooke[_11] * thermal_strain[XX] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + hooke[_11] * thermal_strain[YY]};
            } else if constexpr (std::is_same_v<Expansion, ATE>) {
                return {hooke[_11] * thermal_strain[XX] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + hooke[_11] * thermal_strain[YY],
                    2 * hooke[_66] * thermal_strain[XY]};
            } else static_assert(false, "Unsupported Hooke matrix and thermal expansion tensor combintation.");
        } else if constexpr (std::is_same_v<Hooke, OHM>) {
            using namespace orthotropic_indices;
            if constexpr (std::is_same_v<Expansion, ITE>) {
                return {(hooke[_11] + hooke[_12]) * thermal_strain, 
                        (hooke[_12] + hooke[_22]) * thermal_strain};
            } else if constexpr (std::is_same_v<Expansion, OTE>) {
                return {hooke[_11] * thermal_strain[XX] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + hooke[_22] * thermal_strain[YY]};
            } else if constexpr (std::is_same_v<Expansion, ATE>) {
                return {hooke[_11] * thermal_strain[XX] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + hooke[_22] * thermal_strain[YY],
                    2 * hooke[_66] * thermal_strain[XY]};
            } else static_assert(false, "Unsupported Hooke matrix and thermal expansion tensor combintation.");
        } else if constexpr (std::is_same_v<Hooke, AHM>) {
            using namespace anisotropic_indices;
            if constexpr (std::is_same_v<Expansion, ITE>) {
                return {(hooke[_11] + hooke[_12]) * thermal_strain,
                        (hooke[_12] + hooke[_22]) * thermal_strain,
                        (hooke[_16] + hooke[_26]) * thermal_strain};
            } else if constexpr (std::is_same_v<Expansion, OTE>) {
                return {hooke[_11] * thermal_strain[XX] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + hooke[_22] * thermal_strain[YY],
                        hooke[_16] * thermal_strain[XX] + hooke[_26] * thermal_strain[YY]};
            } else if constexpr (std::is_same_v<Expansion, ATE>) {
                return {hooke[_11] * thermal_strain[XX] + 2 * hooke[_16] * thermal_strain[XY] + hooke[_12] * thermal_strain[YY],
                        hooke[_12] * thermal_strain[XX] + 2 * hooke[_26] * thermal_strain[XY] + hooke[_22] * thermal_strain[YY],
                        hooke[_16] * thermal_strain[XX] + 2 * hooke[_66] * thermal_strain[XY] + hooke[_26] * thermal_strain[YY]};
            } else static_assert(false, "Unsupported Hooke matrix and thermal expansion tensor combintation.");
        } else static_assert(false, "Unsupported Hooke matrix type");
    }

    template<class Hooke, class Expansion>
    std::array<T, 2> operator()(const Hooke& hooke_matrix, const Expansion& thermal_expansion, const size_t e, const size_t i) const {
        using namespace metamath::functions;
        std::array<T, 2> integral = {};
        size_t qshift = _mesh.quad_shift(e);
        const auto& el = _mesh.container().element_2d(e);
        for(const size_t q : el.qnodes()) {
            const auto thermal_stress = calc_stress(hooke_matrix, thermal_expansion, qshift);
            const auto wdN = el.weight(q) * _mesh.derivatives(e, i, q);
            integral[X] += wdN[X] * (thermal_stress[XX] + thermal_stress[XY]);
            integral[Y] += wdN[Y] * (thermal_stress[YY] + thermal_stress[XY]);
            ++qshift;
        }
        return integral;
    }

    template<class Hooke, class Expansion>
    std::array<T, 2> operator()(const Hooke& hooke_matrix, const Expansion& thermal_expansion,
                                const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
                                const size_t eL, const size_t eNL, const size_t iL) const {
        using namespace metamath::functions;
        std::array<T, 2> integral = {};
        const auto& elL = _mesh.container().element_2d(eL);
        const auto& elNL = _mesh.container().element_2d(eNL);
        size_t qshiftL = _mesh.quad_shift(eL);
        for(const size_t qL : elL.qnodes()) {
            std::array<T, 3> inner_integral = {};
            size_t qshiftNL = _mesh.quad_shift(eNL);
            const std::array<T, 2>& qcoordL = _mesh.quad_coord(eL, qL);
            for(const size_t qNL : elNL.qnodes()) {
                const T weight = elNL.weight(qNL) * influence(qcoordL, _mesh.quad_coord(qshiftNL)) * _mesh.jacobian(qshiftNL);
                inner_integral += weight * calc_stress(hooke_matrix, thermal_expansion, qshiftNL);
                ++qshiftNL;
            }
            const auto wdN = elL.weight(qL) * _mesh.derivatives(eL, iL, qL);
            integral[X] += wdN[X] * (inner_integral[XX] + inner_integral[XY]);
            integral[Y] += wdN[Y] * (inner_integral[YY] + inner_integral[XY]);
            ++qshiftL;
        }
        return integral;
    }

public:
    template<class U, class J>
    friend void temperature_condition(Eigen::Matrix<U, Eigen::Dynamic, 1>& f,
                                      const mesh::mesh_2d<U, J>& mesh,
                                      const evaluated_mechanical_parameters<U>& parameters,
                                      const std::vector<U>& delta_temperature);
};

template<class T, class I>
void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                           const mesh::mesh_2d<T, I>& mesh,
                           const evaluated_mechanical_parameters<T>& parameters,
                           const std::vector<T>& delta_temperature) {
    if (delta_temperature.empty())
        return;
    if (delta_temperature.size() != mesh.container().nodes_count())
        throw std::domain_error{"It is impossible to calculate temperature deformations due to "
                                "different dimensions of the mesh and the delta temperature vector."};
    
    const _temperature_condition<T, I> integrator{mesh, delta_temperature};
    const auto process_node = mesh.process_nodes();
#pragma omp parallel for default(none) shared(f, mesh, parameters, integrator, process_node) schedule(dynamic)
    for(size_t node = process_node.front(); node < *process_node.end(); ++node) {
        std::array<T, 2> integral = {};
        for(const I eL : mesh.elements(node)) {
            using namespace metamath::functions;
            const auto& group = mesh.container().group(eL);
            const auto& [model, physical] = parameters.at(group);
            if (physical.thermal_expansion.valueless_by_exception())
                continue;
            const size_t iL = mesh.global_to_local(eL, node);
            if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
                for(const I eNL : mesh.neighbours(eL)) {
                    const auto integrate_nonlocal = [&integrator, &influence = model.influence, eL, eNL, iL]
                                                    (const auto& hooke, const auto& thermal_expansion) {
                        return integrator(hooke, thermal_expansion, influence, eL, eNL, iL);
                    };
                    integral += nonlocal_weight * std::visit(integrate_nonlocal, physical.elastic, physical.thermal_expansion);
                }
            }
            const auto integrate_local = [&integrator, eL, iL](const auto& hooke, const auto& thermal_expansion) {
                return integrator(hooke, thermal_expansion, eL, iL);
            };
            integral += model.local_weight * std::visit(integrate_local, physical.elastic, physical.thermal_expansion);
        }
        f[2 * node + X] += integral[X];
        f[2 * node + Y] += integral[Y];
    }
}

}
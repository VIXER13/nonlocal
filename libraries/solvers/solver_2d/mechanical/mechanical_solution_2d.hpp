#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>
#include <metamath/types/visitor.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T>
class mechanical_solution_2d : public solution_2d<T> {
    using _base = solution_2d<T>;

    std::vector<std::array<T, 2>> _displacement;
    std::vector<std::array<T, 3>> _strain, _stress;

    template<class Hooke_Matrix, class Influence>
    std::array<T, 3> calc_nonlocal_stress(const size_t eL, const Hooke_Matrix& hooke_matrices, 
                                          const std::vector<std::array<T, 3>>& strains, const Influence& influence) const;

    std::vector<std::array<T, 3>> strains_in_quadratures() const;
    void substract_temperature_strains(std::vector<std::array<T, 3>>& strain,
                                       const evaluated_mechanical_parameters<T>& parameters) const;

public:
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh);
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                    const evaluated_mechanical_parameters<T>& parameters, 
                                    std::vector<std::array<T, 2>>&& displacement);
    ~mechanical_solution_2d() noexcept override = default;

    const std::vector<std::array<T, 2>>& displacement() const noexcept;
    const std::vector<std::array<T, 3>>& strain() const noexcept;
    const std::vector<std::array<T, 3>>& stress() const noexcept;

    T calc_energy() const;
    bool is_strain_and_stress_calculated() const noexcept;
    void calc_strain_and_stress(const evaluated_mechanical_parameters<T>& parameters);
};

template<class T>
mechanical_solution_2d<T>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh)
    : _base{mesh}
    , _displacement(mesh->container().nodes_count(), std::array<T, 2>{}) {}

template<class T>
mechanical_solution_2d<T>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                  const evaluated_mechanical_parameters<T>& parameters,
                                                  std::vector<std::array<T, 2>>&& displacement)
    : _base{mesh, {}}
    , _displacement(std::move(displacement)) {}

template<class T>
const std::vector<std::array<T, 2>>& mechanical_solution_2d<T>::displacement() const noexcept {
    return _displacement;
}

template<class T>
const std::vector<std::array<T, 3>>& mechanical_solution_2d<T>::strain() const noexcept {
    return _strain;
}

template<class T>
const std::vector<std::array<T, 3>>& mechanical_solution_2d<T>::stress() const noexcept {
    return _stress;
}

template<class T>
T mechanical_solution_2d<T>::calc_energy() const {
    T integral = 0;
    // if(is_strain_and_stress_calculated()) {
    //     for(size_t e = 0; e < _base::mesh_proxy()->mesh().elements_count(); ++e) {
    //         const auto& el     = _base::mesh_proxy()->mesh().element_2d(e);
    //         const auto J_start = _base::mesh_proxy()->jacobi_matrix(e);
    //         for(size_t i = 0; i < el->nodes_count(); ++i) {
    //             auto J = J_start;
    //             for(size_t q = 0; q < el->qnodes_count(); ++q, ++J) {
    //                 const size_t node = _base::mesh_proxy()->mesh().node_number(e, i);
    //                 integral += el->weight(q) * el->qN(i, q) * mesh::jacobian(*J) *
    //                             (    strain()[XX][node] * stress()[XX][node] +
    //                                  strain()[YY][node] * stress()[YY][node] +
    //                              2 * strain()[XY][node] * stress()[XY][node]);
    //             }
    //         }
    //     }
    // }
    return 0.5 * integral;
}

template<class T>
bool mechanical_solution_2d<T>::is_strain_and_stress_calculated() const noexcept {
    return !strain().empty() && !stress().empty();
}

template<class T>
std::vector<std::array<T, 3>> mechanical_solution_2d<T>::strains_in_quadratures() const {
    const auto& container = _base::mesh().container();
    std::vector<std::array<T, 3>> strains(_base::mesh().quad_shift(container.elements_2d_count()), std::array<T, 3>{});
#pragma omp parallel for
    for(size_t e = 0; e < container.elements_2d_count(); ++e) {
        const auto& el = container.element_2d(e);
        for(size_t q = 0, qshift = _base::mesh().quad_shift(e); q < el.qnodes_count(); ++q, ++qshift) {
            using namespace metamath::operators;
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()}) {
                const auto& disp = displacement()[container.node_number(e, i)];
                const auto& deriv = _base::mesh().derivatives(e, i, q);
                strains[qshift][XX] += disp[X] * deriv[X];
                strains[qshift][YY] += disp[Y] * deriv[Y];
                strains[qshift][XY] += T{0.5} * (disp[X] * deriv[Y] + disp[Y] * deriv[X]);
            }
            strains[qshift] /= _base::mesh().jacobian(qshift);
        }
    }
    return strains;
}

template<class T>
void mechanical_solution_2d<T>::substract_temperature_strains(std::vector<std::array<T, 3>>& strain,
                                                              const evaluated_mechanical_parameters<T>& parameters) const {
    for(const auto& [group, parameter] : parameters)
        std::visit(metamath::types::visitor{
            [](const std::monostate) {},
            [this, &strain, &group](const auto& thermal_strain) {
                using namespace metamath::operators;
                for(const size_t qshift : _base::mesh().quad_shifts(group))
                    strain[qshift] -= thermal_strain[qshift];
            }
        }, parameter.physical.thermal_strain);
}

template<class T>
template<class Hooke_Matrix, class Influence>
std::array<T, 3> mechanical_solution_2d<T>::calc_nonlocal_stress(const size_t eL,
                                                                 const Hooke_Matrix& hooke_matrices,
                                                                 const std::vector<std::array<T, 3>>& strains,
                                                                 const Influence& influence) const {
    std::array<T, 3> nonlocal_stress = {};
    for(const size_t eNL : _base::mesh().neighbours(eL)) {
        const auto& elNL = _base::mesh().container().element_2d(eNL);
        const size_t qshiftNL = _base::mesh().quad_shift(eNL);
        for(const size_t qNL : elNL.qnodes()) {
            const size_t qshift = qshiftNL + qNL;
            const auto& hooke = hooke_matrices.index() ? std::get<Variable>(hooke_matrices)[qshift] : 
                                                         std::get<Constant>(hooke_matrices);
            const T influence_weight = elNL.weight(qNL) * _base::mesh().jacobian(qshift) *
                                       influence(_base::mesh().quad_coord(qshift));
            using namespace metamath::operators;
            nonlocal_stress += influence_weight * calc_stress<T>(hooke, strains[qshift]);
        }
    }
    return nonlocal_stress;
}

template<class T>
void mechanical_solution_2d<T>::calc_strain_and_stress(const evaluated_mechanical_parameters<T>& parameters) {
    if (is_strain_and_stress_calculated())
        return;

    auto strains = strains_in_quadratures();
    _strain = mesh::utils::qnodes_to_nodes(_base::mesh(), strains);
    _stress.resize(strains.size(), std::array<T, 3>{});
    substract_temperature_strains(strains, parameters);
    for(const auto& [group, parameter] : parameters) {
        const auto& [model, physics] = parameter;
        const auto elements = _base::mesh().container().elements(group);
        std::visit([this, &model, &elements, &strains](const auto& hooke_matrices) {
#pragma omp parallel for schedule(dynamic)
            for(size_t eL = elements.front(); eL < *elements.end(); ++eL)
                for(const size_t qshiftL : std::ranges::iota_view{_base::mesh().quad_shift(eL), _base::mesh().quad_shift(eL + 1)}) {
                    using namespace metamath::operators;
                    const auto& hooke = hooke_matrices.index() ? std::get<Variable>(hooke_matrices)[qshiftL] : 
                                                                 std::get<Constant>(hooke_matrices);
                    std::array<T, 3> stress = model.local_weight * calc_stress<T>(hooke, strains[qshiftL]);
                    if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                        const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
                        const auto& qnodeL = _base::mesh().quad_coord(qshiftL);
                        const auto influence = [nonlocal_weight, &influence = model.influence, &qnodeL](const std::array<T, 2>& qnodeNL) {
                            return nonlocal_weight * influence(qnodeL, qnodeNL);
                        };
                        stress += calc_nonlocal_stress(eL, hooke_matrices, strains, influence);
                    }
                    _stress[qshiftL] += stress;
                }
        }, physics.elastic);
    }
    _stress = mesh::utils::qnodes_to_nodes(_base::mesh(), _stress);
    _stress = parallel::all_to_all(_stress, _base::mesh().MPI_ranges());
}

}
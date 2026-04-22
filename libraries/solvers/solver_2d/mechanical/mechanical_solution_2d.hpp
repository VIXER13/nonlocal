#pragma once

#include "mechanical_parameters_2d.hpp"

#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/base/solution_2d.hpp>
#include <metamath/types/visitor.hpp>

namespace nonlocal::solver_2d::mechanical {

template<class T, class I = uint32_t>
class mechanical_solution_2d : public solution_2d<T, I> {
    using _base = solution_2d<T, I>;

    std::array<std::vector<T>, 2> _displacement;
    std::array<std::vector<T>, 3> _strain, _stress;
    std::unordered_map<std::string, evaluated_mechanical_parameters_t<T>> _parameters;

    template<class Hooke_Matrix, class Influence>
    std::array<T, 3> calc_nonlocal_stress(const size_t eL, const Hooke_Matrix& hooke_matrices, 
                                          const std::array<std::vector<T>, 3>& strains, const Influence& influence) const;

    std::array<std::vector<T>, 3> strains_in_quadratures() const;
    void substract_temperature_strains(std::array<std::vector<T>, 3>& strain) const;

public:
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    template<class Vector>
    explicit mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                    const evaluated_mechanical_parameters<T>& parameters, 
                                    const Vector& displacement);
    ~mechanical_solution_2d() noexcept override = default;

    const std::array<std::vector<T>, 2>& displacement() const noexcept;
    const std::array<std::vector<T>, 3>& strain() const noexcept;
    const std::array<std::vector<T>, 3>& stress() const noexcept;

    const evaluated_mechanical_parameters_t<T>& parameters(const std::string& group) const;

    T calc_energy() const;
    bool is_strain_and_stress_calculated() const noexcept;
    void calc_strain_and_stress();
};

template<class T, class I>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(mesh->container().nodes_count(), T{0});
}

template<class T, class I>
template<class Vector>
mechanical_solution_2d<T, I>::mechanical_solution_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                     const evaluated_mechanical_parameters<T>& parameters,
                                                     const Vector& displacement)
    : _base{mesh, get_models(parameters)}
    , _parameters{get_physical_parameters(parameters)} {
    for(std::vector<T>& displacement : _displacement)
        displacement.resize(_base::mesh().container().nodes_count(), T{0});
    for(const size_t i : std::views::iota(0u, _base::mesh().container().nodes_count())) {
        _displacement[X][i] = displacement[2 * i + X];
        _displacement[Y][i] = displacement[2 * i + Y];
    }
}

template<class T, class I>
const std::array<std::vector<T>, 2>& mechanical_solution_2d<T, I>::displacement() const noexcept {
    return _displacement;
}

template<class T, class I>
const std::array<std::vector<T>, 3>& mechanical_solution_2d<T, I>::strain() const noexcept {
    return _strain;
}

template<class T, class I>
const std::array<std::vector<T>, 3>& mechanical_solution_2d<T, I>::stress() const noexcept {
    return _stress;
}

template<class T, class I>
const evaluated_mechanical_parameters_t<T>& mechanical_solution_2d<T, I>::parameters(const std::string& group) const {
    return _parameters.at(group);
}

template<class T, class I>
T mechanical_solution_2d<T, I>::calc_energy() const {
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

template<class T, class I>
bool mechanical_solution_2d<T, I>::is_strain_and_stress_calculated() const noexcept {
    return !strain()[XX].empty() && !strain()[YY].empty() && !strain()[XY].empty() &&
           !stress()[XX].empty() && !stress()[YY].empty() && !stress()[XY].empty();
}

template<class T, class I>
std::array<std::vector<T>, 3> mechanical_solution_2d<T, I>::strains_in_quadratures() const {
    auto [strain11_in_quad, strain12_in_quad] = mesh::utils::gradient_in_qnodes(_base::mesh(), displacement()[X]);
    auto [strain21_in_quad, strain22_in_quad] = mesh::utils::gradient_in_qnodes(_base::mesh(), displacement()[Y]);
    for(const size_t q : std::ranges::iota_view{0zu, strain12_in_quad.size()})
        strain12_in_quad[q] = T{0.5} * (strain12_in_quad[q] + strain21_in_quad[q]);
    return {std::move(strain11_in_quad), std::move(strain22_in_quad), std::move(strain12_in_quad)};
}

template<class T, class I>
void mechanical_solution_2d<T, I>::substract_temperature_strains(std::array<std::vector<T>, 3>& strain) const {
    for(const std::string& group : _base::mesh().container().groups_2d())
        std::visit(metamath::types::visitor{
            [](const std::monostate) {},
            [this, &strain, &group](const auto& thermal_strain) {
                for(const size_t qshift : _base::mesh().quad_shifts(group)) {
                    const auto temperature_strain = thermal_strain[qshift];
                    strain[XX][qshift] -= temperature_strain[XX];
                    strain[YY][qshift] -= temperature_strain[YY];
                    strain[XY][qshift] -= temperature_strain[XY];
                }
            }
        }, parameters(group).thermal_strain);
}

template<class T, class I>
template<class Hooke_Matrix, class Influence>
std::array<T, 3> mechanical_solution_2d<T, I>::calc_nonlocal_stress(const size_t eL,
                                                                    const Hooke_Matrix& hooke_matrices,
                                                                    const std::array<std::vector<T>, 3>& strains,
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
            using namespace metamath::functions;
            nonlocal_stress += influence_weight * calc_stress<T>(hooke, {strains[XX][qshift], strains[YY][qshift], strains[XY][qshift]});
        }
    }
    return nonlocal_stress;
}

template<class T, class I>
void mechanical_solution_2d<T, I>::calc_strain_and_stress() {
    if (is_strain_and_stress_calculated())
        return;

    auto strains = strains_in_quadratures();
    for(const size_t i : std::ranges::iota_view{0zu, 3zu}) {
        _strain[i] = mesh::utils::qnodes_to_nodes(_base::mesh(), strains[i]);
        _stress[i].resize(strains[i].size(), T{0});
    }
    substract_temperature_strains(strains);
    for(const auto& [group, parameter] : _parameters) {
        const auto& model = _base::model(group);
        const auto& hooke_matrices = parameters(group).elastic;
        const auto elements = _base::mesh().container().elements(group);
        std::visit([this, &model, &elements, &strains](const auto& hooke_matrices) {
#pragma omp parallel for schedule(dynamic)
            for(size_t eL = elements.front(); eL < *elements.end(); ++eL)
                for(const size_t qshiftL : std::ranges::iota_view{_base::mesh().quad_shift(eL), _base::mesh().quad_shift(eL + 1)}) {
                    using namespace metamath::functions;
                    const auto& hooke = hooke_matrices.index() ? std::get<Variable>(hooke_matrices)[qshiftL] : 
                                                                 std::get<Constant>(hooke_matrices);
                    std::array<T, 3> stress = model.local_weight * calc_stress<T>(hooke, {strains[XX][qshiftL], strains[YY][qshiftL], strains[XY][qshiftL]});
                    if (theory_type(model.local_weight) == theory_t::NONLOCAL) {
                        const T nonlocal_weight = nonlocal::nonlocal_weight(model.local_weight);
                        const auto& qnodeL = _base::mesh().quad_coord(qshiftL);
                        const auto influence = [nonlocal_weight, &influence = model.influence, &qnodeL](const std::array<T, 2>& qnodeNL) {
                            return nonlocal_weight * influence(qnodeL, qnodeNL);
                        };
                        stress += calc_nonlocal_stress(eL, hooke_matrices, strains, influence);
                    }
                    for(const size_t i : std::ranges::iota_view{0zu, 3zu})
                        _stress[i][qshiftL] += stress[i];
                }
        }, hooke_matrices);

    }
    for(const size_t i : std::ranges::iota_view{0u, 3u}) {
        _stress[i] = mesh::utils::qnodes_to_nodes(_base::mesh(), _stress[i]);
        _stress[i] = parallel::all_to_all(_stress[i], _base::mesh().MPI_ranges());
    }
}

}
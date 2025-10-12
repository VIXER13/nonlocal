#pragma once

#include "evaluate_conductivity.hpp"

#include <solvers/solver_2d/base/matrix_assembler_2d.hpp>

#include <string>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T, std::integral I, std::integral J>
class conductivity_matrix_2d : public matrix_assembler_2d<T, I, J, 1> {
    using _base = matrix_assembler_2d<T, I, J, 1>;

    std::vector<T> _solution; // stub for nonlinear problems

    void create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories,
                                const std::vector<bool>& is_inner, const bool is_symmetric, const bool is_neumann);

    T integrate_basic(const size_t e, const size_t i) const;
    void integral_condition(const bool is_symmetric);

    template<class Material>
    T integrate_local(const Material& conductivity, const size_t e, const size_t i, const size_t j) const;
    template<class Material>
    T integrate_nonlocal(const Material& conductivity, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>)>& influence,
                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

public:
    explicit conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh);
    ~conductivity_matrix_2d() noexcept override = default;

    void compute(const evaluated_conductivity_2d<T>& conductivity, const std::vector<bool>& is_inner,
                 const bool is_symmetric = true, const bool is_neumann = false, const assemble_part part = assemble_part::FULL);
};

template<std::floating_point T, std::integral I, std::integral J>
conductivity_matrix_2d<T, I, J>::conductivity_matrix_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh)
    : _base{mesh} {}

template<std::floating_point T, std::integral I, std::integral J>
void conductivity_matrix_2d<T, I, J>::create_matrix_portrait(const std::unordered_map<std::string, theory_t> theories, 
                                                             const std::vector<bool>& is_inner,
                                                             const bool is_symmetric,
                                                             const bool is_neumann) {
    const size_t cols = _base::cols() + is_neumann;
    const size_t rows = _base::rows() == _base::cols() ?
                        cols : _base::rows() + (is_neumann && parallel::is_last_process());
    _base::matrix().inner().resize(rows, cols);
    _base::matrix().bound().resize(rows, cols);
    if (is_neumann) {
        for(const size_t row : std::views::iota(0u, rows))
            _base::matrix().inner().outerIndexPtr()[row + 1] = 1;
        if (!is_symmetric && parallel::is_last_process())
            _base::matrix().inner().outerIndexPtr()[rows] = cols - 1;
    }
    _base::init_shifts(theories, is_inner, is_symmetric);
    static constexpr bool SORT_INDICES = false;
    _base::init_indices(theories, is_inner, is_symmetric, SORT_INDICES);
    if (is_neumann) {
        for(const size_t row : std::ranges::iota_view{0u, rows}) {
            const size_t index = _base::matrix().inner().outerIndexPtr()[row + 1] - 1;
            _base::matrix().inner().innerIndexPtr()[index] = _base::mesh().container().nodes_count();
        }
        if (!is_symmetric && parallel::is_last_process()) 
            for(const size_t col : std::ranges::iota_view{0u, cols}) {
                const size_t index = _base::matrix().inner().outerIndexPtr()[rows - 1] + col;
                _base::matrix().inner().innerIndexPtr()[index] = col;
            }
    }
    utils::sort_indices(_base::matrix().inner());
    utils::sort_indices(_base::matrix().bound());
    logger::info() << "Matrix portrait is formed" << std::endl;
}

template<std::floating_point T, std::integral I, std::integral J>
T conductivity_matrix_2d<T, I, J>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes())
        integral += el.weight(q) * el.qN(i, q) * _base::mesh().jacobian(e, q);
    return integral;
}

template<std::floating_point T, std::integral I, std::integral J>
void conductivity_matrix_2d<T, I, J>::integral_condition(const bool is_symmetric) {
    const auto process_nodes = _base::rows() == _base::cols() ?
                               std::ranges::iota_view<size_t, size_t>{0u, size_t(_base::matrix().inner().cols()) - 1} :
                               std::get<std::ranges::iota_view<size_t, size_t>>(_base::nodes_for_processing());
#pragma omp parallel for default(none) shared(process_nodes, is_symmetric)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix().inner().coeffRef(node - process_nodes.front(), _base::mesh().container().nodes_count());
        for(const I e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        if (!is_symmetric && parallel::is_last_process())
            _base::matrix().inner().coeffRef(_base::matrix().inner().rows() - 1, node) = val;
    }
    if (!is_symmetric && parallel::MPI_size() > 1 && parallel::is_last_process()) {
#pragma omp parallel for default(none) shared(process_nodes, is_symmetric)
        for(size_t node = 0; node < process_nodes.front(); ++node) {
            T& val = _base::matrix().inner().coeffRef(_base::matrix().inner().rows() - 1, node);
            for(const I e : _base::mesh().elements(node))
                val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        }
    }
}

template<std::floating_point T, std::integral I, std::integral J>
template<class Material>
T conductivity_matrix_2d<T, I, J>::integrate_local(const Material& conductivity, const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const size_t qshift = _base::mesh().quad_shift(e);
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes()) {
        const auto& dNi = _base::mesh().derivatives(e, i, q);
        const auto& dNj = _base::mesh().derivatives(e, j, q);
        T value = T{0};
        const auto& conduct = conductivity.index() ? std::get<1>(conductivity)[qshift + q] : std::get<0>(conductivity);
        if constexpr (std::is_same_v<Material, evaluated_isotropic_conductivity_t<T>>)
            value = conduct * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]);
        else if constexpr (std::is_same_v<Material, evaluated_orthotropic_conductivity_t<T>>)
            value = conduct[X] * dNi[X] * dNj[X] + conduct[Y] * dNi[Y] * dNj[Y];
        else if constexpr (std::is_same_v<Material, evaluated_anisotropic_conductivity_t<T>>)
            value = dNi[X] * (conduct[X][X] * dNj[X] + conduct[X][Y] * dNj[Y]) + 
                    dNi[Y] * (conduct[Y][X] * dNj[X] + conduct[Y][Y] * dNj[Y]);
        integral += el.weight(q) * value / _base::mesh().jacobian(e, q);
    }
    return integral;
}

template<std::floating_point T, std::integral I, std::integral J>
template<class Material>
T conductivity_matrix_2d<T, I, J>::integrate_nonlocal(const Material& conductivity, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>)>& influence,
                                                      const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const {
    T integral = T{0};
    const auto& elL  = _base::mesh().container().element_2d(eL );
    const auto& elNL = _base::mesh().container().element_2d(eNL);
    const size_t qshiftNL = _base::mesh().quad_shift(eNL);
    for(const size_t qL : elL.qnodes()) {
        std::array<T, 2> inner_integral = {};
        const auto& qcoordL = _base::mesh().quad_coord(eL, qL);
        const auto& dNi = _base::mesh().derivatives(eL, iL, qL);
        for(const size_t qNL : elNL.qnodes()) {
            const auto& conduct = conductivity.index() ? std::get<1>(conductivity)[qshiftNL + qNL] : std::get<0>(conductivity);
            const auto& dNj = _base::mesh().derivatives(eNL, jNL, qNL);
            const T influence_weight = elNL.weight(qNL) * influence(qcoordL, _base::mesh().quad_coord(eNL, qNL));
            if constexpr (std::is_same_v<Material, isotropic_conductivity_t<T>>) {
                using namespace metamath::functions;
                inner_integral += influence_weight * conduct * dNj;
            } else if constexpr (std::is_same_v<Material, orthotropic_conductivity_t<T>>) {
                inner_integral[X] += influence_weight * conduct[X] * dNj[X];
                inner_integral[Y] += influence_weight * conduct[Y] * dNj[Y];
            } else if constexpr (std::is_same_v<Material, anisotropic_conductivity_t<T>>) {
                inner_integral[X] += influence_weight * (conduct[X][X] * dNj[X] + conduct[X][Y] * dNj[Y]);
                inner_integral[Y] += influence_weight * (conduct[Y][X] * dNj[X] + conduct[Y][Y] * dNj[Y]);
            }
        }
        integral += elL.weight(qL) * (dNi[X] * inner_integral[X] + dNi[Y] * inner_integral[Y]);
    }
    return integral;
}

template<std::floating_point T, std::integral I, std::integral J>
void conductivity_matrix_2d<T, I, J>::compute(const evaluated_conductivity_2d<T>& conductivity, const std::vector<bool>& is_inner,
                                              const bool is_symmetric, const bool is_neumann, const assemble_part part) {
    logger::info() << "Thermal conductivity matrix assembly started" << std::endl;
    const std::unordered_map<std::string, theory_t> theories = part == assemble_part::LOCAL ? 
                                                                local_theories(_base::mesh().container()) : 
                                                                theories_types(conductivity);
    create_matrix_portrait(theories, is_inner, is_symmetric, is_neumann);
    if (is_neumann)
        integral_condition(is_symmetric);
    _base::calc_coeffs(theories, is_inner, is_symmetric,
        [this, &conductivity](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& group_params = conductivity.at(group);
            const auto& model = group_params.model;
            const auto& physic = group_params.physical;
            const T integral = std::visit(metamath::visitor{
                [&](const evaluated_isotropic_conductivity_t<T>& conductivity) { return integrate_local(conductivity, e, i, j); },
                [&](const evaluated_orthotropic_conductivity_t<T>& conductivity) { return integrate_local(conductivity, e, i, j); },
                [&](const evaluated_anisotropic_conductivity_t<T>& conductivity) { return integrate_local(conductivity, e, i, j); }
            }, physic);
            return model.local_weight * integral;
        },
        [this, &conductivity](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& group_params = conductivity.at(group);
            const auto& model = group_params.model;
            const auto& physic = group_params.physical;
            const T integral = std::visit(metamath::visitor{
                [&](const evaluated_isotropic_conductivity_t<T>& conductivity) { return integrate_nonlocal(conductivity, model.influence, eL, eNL, iL, jNL); },
                [&](const evaluated_orthotropic_conductivity_t<T>& conductivity) { return integrate_nonlocal(conductivity, model.influence, eL, eNL, iL, jNL); },
                [&](const evaluated_anisotropic_conductivity_t<T>& conductivity) { return integrate_nonlocal(conductivity, model.influence, eL, eNL, iL, jNL); }
            }, physic);
            return nonlocal::nonlocal_weight(model.local_weight) * integral;
        }
    );
    logger::info() << "Thermal conductivity matrix assembly finished" << std::endl;
}

}
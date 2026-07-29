#pragma once

#include "evaluate_conductivity.hpp"

#include <solvers/solver_2d/base/matrix_assembler.hpp>

#include <string>
#include <iostream>

namespace nonlocal::solver_2d::thermal {
    
template<std::floating_point T>
class conductivity_matrix_2d : public matrix_assembler_base<T> {
    using _base = matrix_assembler_base<T>;

    void create_matrix_portrait(const problem_settings& settings);

    T integrate_basic(const size_t e, const size_t i) const;
    void integral_condition(const bool is_symmetric);

    template<class Conductivity>
    T integrate_local(const Conductivity& conductivity, const size_t e, const size_t i, const size_t j) const;
    template<class Conductivity>
    T integrate_nonlocal(const Conductivity& conductivity, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) const;

public:
    explicit conductivity_matrix_2d(const mesh::mesh_2d<T>& mesh);
    ~conductivity_matrix_2d() noexcept override = default;

    void compute(const evaluated_conductivity_2d<T>& conductivity, const problem_settings& settings);
};

template<std::floating_point T>
conductivity_matrix_2d<T>::conductivity_matrix_2d(const mesh::mesh_2d<T>& mesh)
    : _base{mesh} {}

template<std::floating_point T>
void conductivity_matrix_2d<T>::create_matrix_portrait(const problem_settings& settings) {
    const size_t cols = _base::mesh().container().nodes_count() + settings.is_neumann;
    const size_t rows = _base::rows() + (settings.is_neumann && parallel::is_last_process());
    _base::matrix().portrait.set_size(rows, cols);

    if (settings.is_neumann) {
        for(const size_t row : std::views::iota(0u, rows))
            _base::matrix().portrait.shifts[row + 1] = 1;
        if (!settings.is_symmetric() && parallel::is_last_process())
            _base::matrix().portrait.shifts[rows] = cols - 1;
    }

    _base::init_shifts(settings);
    static constexpr bool Sort_Indices = false;
    _base::init_indices(settings, Sort_Indices);

    if (settings.is_neumann) {
        for(const size_t row : std::ranges::iota_view{0u, rows}) {
            const size_t index = _base::matrix().portrait.shifts[row + 1] - 1;
            _base::matrix().portrait.indices[index] = _base::mesh().container().nodes_count();
        }
        if (!settings.is_symmetric() && parallel::is_last_process())
            for(const size_t col : std::ranges::iota_view{0u, cols}) {
                const size_t index = _base::matrix().portrait.indices[rows - 1] + col;
                _base::matrix().portrait.indices[index] = col;
            }
    }
    
    _base::matrix().portrait.sort_indices();
    logger::info() << "Matrix portrait is formed" << std::endl;
}

template<std::floating_point T>
T conductivity_matrix_2d<T>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes())
        integral += el.weight(q) * el.qN(i, q) * _base::mesh().jacobian(e, q);
    return integral;
}

template<std::floating_point T>
void conductivity_matrix_2d<T>::integral_condition(const bool is_symmetric) {
    const auto process_nodes = std::get<std::ranges::iota_view<size_t, size_t>>(_base::processing_nodes());

#pragma omp parallel for
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        T& val = _base::matrix()(node - process_nodes.front(), _base::mesh().container().nodes_count());
        for(const size_t e : _base::mesh().elements(node))
            val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        if (!is_symmetric && parallel::is_last_process())
            _base::matrix()(_base::matrix().rows() - 1, node) = val;
    }

    if (!is_symmetric && parallel::MPI_size() > 1 && parallel::is_last_process()) {
#pragma omp parallel for
        for(size_t node = 0; node < process_nodes.front(); ++node) {
            T& val = _base::matrix()(_base::matrix().rows() - 1, node);
            for(const size_t e : _base::mesh().elements(node))
                val += integrate_basic(e, _base::mesh().global_to_local(e, node));
        }
    }
}

template<std::floating_point T>
template<class Conductivity>
T conductivity_matrix_2d<T>::integrate_local(const Conductivity& conductivity, const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const size_t qshift = _base::mesh().quad_shift(e);
    const auto& el = _base::mesh().container().element_2d(e);
    for(const size_t q : el.qnodes()) {
        const auto& dNi = _base::mesh().derivatives(e, i, q);
        const auto& dNj = _base::mesh().derivatives(e, j, q);
        T value = T{0};
        const auto& conduct = conductivity.index() ? std::get<Variable>(conductivity)[qshift + q] :
                                                     std::get<Constant>(conductivity);
        if constexpr (std::is_same_v<Conductivity, evaluated_isotropic_conductivity_t<T>>)
            value = conduct * (dNi[X] * dNj[X] + dNi[Y] * dNj[Y]);
        else if constexpr (std::is_same_v<Conductivity, evaluated_orthotropic_conductivity_t<T>>)
            value = conduct[X] * dNi[X] * dNj[X] + conduct[Y] * dNi[Y] * dNj[Y];
        else if constexpr (std::is_same_v<Conductivity, evaluated_anisotropic_conductivity_t<T>>)
            value = dNi[X] * (conduct[XX] * dNj[X] + conduct[XY] * dNj[Y]) +
                    dNi[Y] * (conduct[XY] * dNj[X] + conduct[YY] * dNj[Y]);
        else
            static_assert(false, "Unsupported coefficient type.");
        integral += el.weight(q) * value / _base::mesh().jacobian(e, q);
    }
    return integral;
}

template<std::floating_point T>
template<class Conductivity>
T conductivity_matrix_2d<T>::integrate_nonlocal(const Conductivity& conductivity, const std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)>& influence,
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
            const auto& conduct = conductivity.index() ? std::get<Variable>(conductivity)[qshiftNL + qNL] :
                                                         std::get<Constant>(conductivity);
            const auto& dNj = _base::mesh().derivatives(eNL, jNL, qNL);
            const T influence_weight = elNL.weight(qNL) * influence(qcoordL, _base::mesh().quad_coord(eNL, qNL));
            if constexpr (std::is_same_v<Conductivity, evaluated_isotropic_conductivity_t<T>>) {
                using namespace metamath::operators;
                inner_integral += influence_weight * conduct * dNj;
            } else if constexpr (std::is_same_v<Conductivity, evaluated_orthotropic_conductivity_t<T>>) {
                inner_integral[X] += influence_weight * conduct[X] * dNj[X];
                inner_integral[Y] += influence_weight * conduct[Y] * dNj[Y];
            } else if constexpr (std::is_same_v<Conductivity, evaluated_anisotropic_conductivity_t<T>>) {
                inner_integral[X] += influence_weight * (conduct[XX] * dNj[X] + conduct[XY] * dNj[Y]);
                inner_integral[Y] += influence_weight * (conduct[XY] * dNj[X] + conduct[YY] * dNj[Y]);
            } else
                static_assert(false, "Unsupported coefficient type.");
        }
        integral += elL.weight(qL) * (dNi[X] * inner_integral[X] + dNi[Y] * inner_integral[Y]);
    }
    return integral;
}

template<std::floating_point T>
void conductivity_matrix_2d<T>::compute(const evaluated_conductivity_2d<T>& conductivity, const problem_settings& settings) {
    logger::info() << "Thermal conductivity matrix assembly started" << std::endl;
    _base::matrix().clear();
    create_matrix_portrait(settings);
    if (settings.is_neumann)
        integral_condition(settings.is_symmetric());
    _base::calc_coeffs(settings,
        [this, &conductivity](const std::string& group, const size_t e, const size_t i, const size_t j) {
            const auto& [model, physic] = conductivity.at(group);
            return model.local_weight * std::visit([this, e, i, j](const auto& conductivity) {
                return integrate_local(conductivity, e, i, j);
            }, physic);
        },
        [this, &conductivity](const std::string& group, const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
            const auto& [model, physic] = conductivity.at(group);
            return nonlocal::nonlocal_weight(model.local_weight) * std::visit([this, &model, eL, eNL, iL, jNL](const auto& conductivity) {
                return integrate_nonlocal(conductivity, model.influence, eL, eNL, iL, jNL);
            }, physic);
        }
    );
    logger::info() << "Thermal conductivity matrix assembly finished" << std::endl;
}

}
#pragma once

#include "thermal_boundary_conditions_1d.hpp"

#include <solvers/solver_1d/base/boundary_condition_nonlinear_1d.hpp>

#include <ranges>

namespace nonlocal::solver_1d::thermal {

template<bool Is_Stationary, class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part,
                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev, const T time_step = T{ 1 }) {
    using metamath::constants::Stefan_Boltzmann_Constant;
    using metamath::functions::power;
    const std::array<size_t, 2> indices = { 0, size_t(matrix.rows() - 1) };
    for (const size_t b : std::ranges::iota_view{ 0zu, 2zu }) {
        if (const auto* const condition = dynamic_cast<const radiation_1d<T>*>(boundaries_conditions[b].get())) {
            const size_t index = indices[b];
            const T previous = temperature_prev[index];
            const T value = condition->emissivity() * Stefan_Boltzmann_Constant<T> * power<4>(previous);
            const T derivative = T{ 4 } * Stefan_Boltzmann_Constant<T> * condition->emissivity() * power<3>(previous);
            boundary_condition_nonlinear_1d<Is_Stationary>(matrix, right_part, index, previous, value, derivative, time_step);
        }
    }
}

} // namespace nonlocal::solver_1d::thermal

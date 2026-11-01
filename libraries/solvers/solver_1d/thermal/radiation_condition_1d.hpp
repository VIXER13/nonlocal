#pragma once

#include "thermal_boundary_conditions_1d.hpp"

#include <solvers/solver_1d/base/boundary_condition_nonlinear_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<bool Is_Stationary, class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part,
                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev, const T time_step = T{ 1 }) {
    boundary_condition_nonlinear_1d<Is_Stationary>(
        matrix, right_part, normalize_boundaries_1d(boundaries_conditions), temperature_prev, time_step);
}

} // namespace nonlocal::solver_1d::thermal

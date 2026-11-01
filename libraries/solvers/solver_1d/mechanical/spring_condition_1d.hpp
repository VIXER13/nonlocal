#pragma once

#include "mechanical_boundary_conditions_1d.hpp"

#include <solvers/solver_1d/base/boundary_condition_robin_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<class T, class I>
void spring_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                             const mechanical_boundary_condition_1d<T>& boundary_condition, const size_t index) {
    boundary_condition_robin_1d(matrix, boundary_condition.data(), index);
}

template<class T, class I>
void spring_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                             const mechanical_boundaries_conditions_1d<T>& boundaries_conditions) {
    boundary_condition_robin_1d(matrix, normalize_boundaries_1d(boundaries_conditions));
}

} // namespace nonlocal::solver_1d::mechanical

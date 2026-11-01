#pragma once

#include "boundary_data_1d.hpp"

#include <Eigen/Sparse>

namespace nonlocal::solver_1d {

template<class T, class I>
void boundary_condition_robin_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                 const boundary_data_1d<T>& condition, const size_t index) {
    if (condition.robin_coefficient)
        matrix.coeffRef(index, index) += *condition.robin_coefficient;
}

template<class T, class I>
void boundary_condition_robin_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                 const boundaries_data_1d<T>& conditions) {
    boundary_condition_robin_1d(matrix, conditions.front(), 0);
    boundary_condition_robin_1d(matrix, conditions.back(), size_t(matrix.rows() - 1));
}

} // namespace nonlocal::solver_1d

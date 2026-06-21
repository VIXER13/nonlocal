#pragma once

#include "mechanical_boundary_conditions_1d.hpp"
#include <Eigen/Sparse>
#include <ranges>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T, std::integral I>
void spring_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                         const mechanical_boundary_condition_1d<T>& boundary_condition,
                         const size_t index) {
    if (const auto* const condition = dynamic_cast<const spring_1d<T>*>(&boundary_condition))
        matrix.coeffRef(index, index) += condition->stiffness();
}

template<std::floating_point T, std::integral I>
void spring_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                         const mechanical_boundaries_conditions_1d<T>& boundaries_conditions) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        spring_condition_1d(matrix, *boundaries_conditions[b], indices[b]);
}

}
#pragma once

#include "thermal_boundary_conditions_1d.hpp"
#include <Eigen/Sparse>
#include <ranges>

namespace nonlocal::solver_1d::thermal {

template<class T, class I>
void convection_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                             const thermal_boundary_condition_1d<T>& boundary_condition,
                             const size_t index) {
    if (const auto* const condition = dynamic_cast<const convection_1d<T>*>(&boundary_condition))
        matrix.coeffRef(index, index) += condition->heat_transfer();
}

template<class T, class I>
void convection_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                             const thermal_boundaries_conditions_1d<T>& boundaries_conditions) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        convection_condition_1d(matrix, *boundaries_conditions[b], indices[b]);
}

}
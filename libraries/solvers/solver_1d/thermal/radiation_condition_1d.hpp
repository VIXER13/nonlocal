#ifndef NONLOCAL_RADIATION_CONDITION_1D_HPP
#define NONLOCAL_RADIATION_CONDITION_1D_HPP

#include "thermal_boundary_conditions_1d.hpp"
#include <eigen3/Eigen/Sparse>
#include <ranges>

namespace nonlocal::thermal {

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            const thermal_boundary_condition_1d<T>& boundary_condition,
                            const size_t index) {
    if (const auto* const condition = dynamic_cast<const convection_1d<T>*>(&boundary_condition))
        matrix.coeffRef(index, index) += condition -> matrix_value();
}

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                            const boundaries_conditions_1d<T>& boundaries_conditions) {
    const std::array<size_t, 2> indices = {0, size_t(matrix.outerSize() - 1)};
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        radiation_condition_1d(matrix, *boundaries_conditions[b], indices[b]);
}

}

#endif
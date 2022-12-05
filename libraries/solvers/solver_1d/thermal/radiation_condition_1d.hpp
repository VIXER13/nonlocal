#ifndef NONLOCAL_RADIATION_CONDITION_1D_HPP
#define NONLOCAL_RADIATION_CONDITION_1D_HPP

#include "thermal_boundary_conditions_1d.hpp"
#include <eigen3/Eigen/Sparse>
#include <ranges>

namespace nonlocal::thermal {

template<class T, class I>
void radiation_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                             const thermal_boundary_condition_1d& boundary_condition,
                             const size_t index) {
    if (const auto* const condition = dynamic_cast<const convection_1d<T>*>(&boundary_condition))
        matrix.coeffRef(index, index) += condition -> matrix_value();
}

}

#endif
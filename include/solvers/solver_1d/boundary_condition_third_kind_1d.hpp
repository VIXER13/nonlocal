#ifndef BOUNDARY_CONDITION_THIRD_KIND_1D_HPP
#define BOUNDARY_CONDITION_THIRD_KIND_1D_HPP

#include "../solvers_constants.hpp"
#include "boundary_condition_1d.hpp"
#include <eigen3/Eigen/Sparse>

namespace nonlocal {

template<class T, class I>
void boundary_condition_third_kind_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                      const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                      const std::array<T, 2>& alpha) {
    if (boundary_condition.front().first == boundary_condition_t::THIRD_KIND)
        K.coeffRef(0, 0) -= alpha.front();
    if (boundary_condition.back().first  == boundary_condition_t::THIRD_KIND)
        K.coeffRef(K.rows()-1, K.cols()-1) -= alpha.back();
}

}

#endif
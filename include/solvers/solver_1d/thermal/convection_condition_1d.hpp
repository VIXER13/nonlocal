#ifndef CONVECTION_CONDITION_1D_HPP
#define CONVECTION_CONDITION_1D_HPP

#include "../../solvers_constants.hpp"
#include "boundary_condition_1d.hpp"
#include <eigen3/Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I>
void convection_condition_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                             const std::array<boundary_condition_t, 2> bound_cond,
                             const std::array<T, 2>& alpha) {
    if (bound_cond.front() == boundary_condition_t::CONVECTION)
        K.coeffRef(0, 0) -= alpha.front();
    if (bound_cond.back()  == boundary_condition_t::CONVECTION)
        K.coeffRef(K.rows()-1, K.cols()-1) -= alpha.back();
}

}

#endif
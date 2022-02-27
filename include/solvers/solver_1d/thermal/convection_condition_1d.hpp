#ifndef CONVECTION_CONDITION_1D_HPP
#define CONVECTION_CONDITION_1D_HPP

#include "../../solvers_constants.hpp"
#include "boundary_condition_1d.hpp"
#include <eigen3/Eigen/Sparse>
#include <optional>

namespace nonlocal::thermal {

template<class T, class I>
void convection_condition_matrix_part_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                         const std::array<boundary_condition_t, 2> boundary_condition,
                                         const std::array<T, 2>& alpha) {
    if (boundary_condition.front() == boundary_condition_t::CONVECTION)
        K.coeffRef(0, 0) -= alpha.front();
    if (boundary_condition.back()  == boundary_condition_t::CONVECTION)
        K.coeffRef(K.rows()-1, K.cols()-1) -= alpha.back();
}

template<class T>
void convection_condition_right_part_1d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                        const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                        const std::array<T, 2>& alpha) {
    if (boundary_condition.front().type == boundary_condition_t::CONVECTION)
        f[0] += alpha.front() * boundary_condition.front().val;
    if (boundary_condition.back().type  == boundary_condition_t::CONVECTION)
        f[f.size()-1] += alpha.back() * boundary_condition.back().val;
}

}

#endif
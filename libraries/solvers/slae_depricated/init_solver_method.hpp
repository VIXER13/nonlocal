#pragma once

#include "conjugate_gradient.hpp"
#include "eigen_preconditioner.hpp"
#include "stable_biconjugate_gradient.hpp"

#include <concepts>
#include <memory>

namespace nonlocal::slae {

template<std::floating_point T, std::integral I>
std::unique_ptr<iterative_solver_base<T, I>> init_iterative_solver(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, const bool is_symmetric) {
    if (is_symmetric)
        return std::make_unique<conjugate_gradient<T, I>>(matrix);
    return std::make_unique<stable_biconjugate_gradient<T, I>>(matrix);
}

template<std::floating_point T, std::integral I>
std::unique_ptr<preconditioner_base<T, I>> init_preconditioner(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, const bool is_symmetric) {
    std::unique_ptr<preconditioner_base<T, I>> preconditioner;
    if (is_symmetric)
        preconditioner = std::make_unique<slae::eigen_ILLT_preconditioner<T, I>>(matrix);
    else
        preconditioner = std::make_unique<slae::eigen_ILUT_preconditioner<T, I>>(matrix);
    if (preconditioner->computation_info() == Eigen::Success)
        return preconditioner;
    return nullptr;
}

}
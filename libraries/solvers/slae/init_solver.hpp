#pragma once

#include "conjugate_gradient.hpp"
#include "ildlt_preconditioner.hpp"
#include "ilu0_preconditioner.hpp"
#include "stable_biconjugate_gradient.hpp"

#include <memory>

namespace nonlocal::slae {

template<class T, std::integral I, std::integral J>
std::unique_ptr<iterative_solver_base<T, I, J>> init_iterative_solver(const metamath::linear::sparse_matrix<T, I, J>& matrix, const bool is_symmetric) {
    if (is_symmetric)
        return std::make_unique<conjugate_gradient<T, I, J>>(matrix);
    return std::make_unique<stable_biconjugate_gradient<T, I, J>>(matrix);
}

template<std::floating_point T, std::integral I, std::integral J>
std::unique_ptr<preconditioner_base<T>> init_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix, const bool is_symmetric) {
    if (is_symmetric)
        return std::make_unique<ildlt_preconditioner<T, I, J>>(std::move(matrix));
    return std::make_unique<ilu0_preconditioner<T, I, J>>(std::move(matrix));
}

}
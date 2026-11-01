#pragma once

#include "boundary_data_1d.hpp"

#include <Eigen/Sparse>

#include <ranges>

namespace nonlocal::solver_1d {

// Linearize the boundary law at v: g(u) ~= g(v) + g'(v)*(u - v).
// A and f below denote the system before this boundary contribution; products act on boundary nodes.
// Stationary Newton step (v is the current iterate, dt defaults to 1):
//   R = A*v - f - dt*g(v), J = A - dt*diag(g'(v)); the caller updates u = v - J^(-1)*R.
// Transient step (v is the previous time layer):
//   (A + dt*diag(g'(v)))*u = f + dt*(g'(v)*v - g(v)).
// Thus g enters the stationary residual with a minus sign and the transient equation with a plus sign.
// These signs preserve the existing solver conventions; callbacks return g and g' themselves.
template<bool Is_Stationary, class T, class I>
void boundary_condition_nonlinear_1d(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix,
                                     Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part,
                                     const boundaries_data_1d<T>& conditions,
                                     const Eigen::Matrix<T, Eigen::Dynamic, 1>& solution_prev,
                                     const T time_step = T{ 1 }) {
    const std::array<size_t, 2> indices = { 0, size_t(matrix.rows() - 1) };
    for (const size_t b : std::ranges::iota_view{ 0zu, 2zu }) {
        if (!conditions[b].nonlinear)
            continue;
        const size_t index = indices[b];
        const auto& op = *conditions[b].nonlinear;
        const T previous = solution_prev[index];
        const T value = op.value(previous);
        const T derivative = op.derivative(previous);
        if constexpr (Is_Stationary) {
            // The stationary solver passes its residual, rather than the load vector.
            matrix.coeffRef(index, index) -= time_step * derivative;
            right_part[index] -= time_step * value;
        } else {
            matrix.coeffRef(index, index) += time_step * derivative;
            right_part[index] += time_step * (derivative * previous - value);
        }
    }
}

} // namespace nonlocal::solver_1d

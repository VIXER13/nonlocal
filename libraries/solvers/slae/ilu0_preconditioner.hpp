#pragma once

#include "preconditioner_base.hpp"

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/types/traits.hpp>
#include <metamath/utils/operators.hpp>

#include <vector>
#include <iostream>

namespace nonlocal::slae {

// ILU0 (Incomplete LU with zero fill-in) preconditioner for general sparse matrices.
// T may be scalar or square_matrix<floating_point, N> for block structure.
// Factorization: A ≈ L U, sparsity pattern preserved.
// L is stored in the lower triangle (unit diagonal implicit), U in the upper triangle including the diagonal.
template<class T, std::integral I, std::integral J>
class ilu0_preconditioner final : public preconditioner_base<T> {
    using typename preconditioner_base<T>::entity_t;

    metamath::linear::sparse_matrix<T, I, J> _matrix;

    template<std::floating_point U>
    static void compute(metamath::linear::sparse_matrix<U, I, J>& matrix) {
        for(const size_t i : std::ranges::iota_view{0zu, matrix.rows()}) {
            // Process lower entries of row i: for each k < i with (i,k) in pattern
            for(const size_t s : matrix.portrait.shifts_range(i))
                if (const size_t k = matrix.portrait.indices[s]; k < i) {
                    using namespace metamath::linear;
                    // l[i,k] = a[i,k] * inv(u[k,k])
                    matrix.values[s] *= matrix(k, k);
                    const auto& lik = matrix.values[s];
                    // Update all remaining entries in row i after position k
                    for(const size_t si : matrix.portrait.shifts_range(i))
                        if (const size_t j = matrix.portrait.indices[si]; j > k && matrix.portrait.contains(k, j))
                            matrix.values[si] -= lik * matrix(k, j);
                }
            matrix(i, i) = metamath::linear::inverse(matrix(i, i));
        }
    }

    template<std::floating_point U, size_t N>
    static void compute(metamath::linear::sparse_matrix<metamath::linear::square_matrix<U, N>, I, J>& matrix) {
        for(const size_t i : std::ranges::iota_view{0zu, matrix.rows()})
            for(const size_t i_dof : std::ranges::iota_view{0zu, N}) {
                for(const size_t s : matrix.portrait.shifts_range(i))
                    for(const size_t k_dof : std::ranges::iota_view{0zu, N})
                        if (const size_t k = matrix.portrait.indices[s]; N * k + k_dof < N * i + i_dof) {
                            using namespace metamath::linear;
                            matrix.values[s][i_dof][k_dof] *= matrix(k, k)[k_dof][k_dof];
                            const auto& lik = matrix.values[s][i_dof][k_dof];
                            for(const size_t si : matrix.portrait.shifts_range(i))
                                if (const size_t j = matrix.portrait.indices[si]; matrix.portrait.contains(k, j))
                                    for(const size_t j_dof : std::ranges::iota_view{0zu, N})
                                        if (N * j + j_dof > N * k + k_dof)
                                            matrix.values[si][i_dof][j_dof] -= lik * matrix(k, j)[k_dof][j_dof];
                        }
                matrix(i, i)[i_dof][i_dof] = metamath::linear::inverse(matrix(i, i)[i_dof][i_dof]);
            }
    }

    template<std::floating_point U>
    static std::vector<entity_t> solve(const metamath::linear::sparse_matrix<U, I, J>& matrix, const std::vector<entity_t>& rhs) {
        std::vector<entity_t> result = rhs;

        // Forward substitution: L z = rhs  (L has implicit unit diagonal)
        // z[i] = rhs[i] - sum_{k < i, (i,k) in pattern} l[i,k] * z[k]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            for(const size_t s : matrix.portrait.shifts_range(i))
                if (const size_t k = matrix.portrait.indices[s]; k < i)
                    result[i] -= matrix.values[s] * result[k];

        // Backward substitution: U x = z
        // x[i] = inv(u[i,i]) * (z[i] - sum_{j > i, (i,j) in pattern} u[i,j] * x[j])
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()} | std::views::reverse) {
            for(const size_t s : matrix.portrait.shifts_range(i))
                if (const size_t j = matrix.portrait.indices[s]; j > i)
                    result[i] -= matrix.values[s] * result[j];
            result[i] = matrix(i, i) * result[i];
        }

        return result;
    }

    template<std::floating_point U, size_t N>
    static std::vector<entity_t> solve(const metamath::linear::sparse_matrix<metamath::linear::square_matrix<U, N>, I, J>& matrix,
                                       const std::vector<entity_t>& rhs) {
        std::vector<entity_t> result = rhs;

        // Forward substitution: L z = rhs  (L has implicit unit diagonal)
        // z[i] = rhs[i] - sum_{k < i, (i,k) in pattern} l[i,k] * z[k]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            for(const size_t i_dof : std::ranges::iota_view{0zu, N})
                for(const size_t s : matrix.portrait.shifts_range(i)) 
                    for(const size_t k_dof : std::ranges::iota_view{0zu, N})
                        if (const size_t k = matrix.portrait.indices[s]; N * k + k_dof < N * i + i_dof)
                            result[i][i_dof] -= matrix.values[s][i_dof][k_dof] * result[k][k_dof];

        // Backward substitution: U x = z
        // x[i] = inv(u[i,i]) * (z[i] - sum_{j > i, (i,j) in pattern} u[i,j] * x[j])
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()} | std::views::reverse) {
            for(const size_t i_dof : std::ranges::iota_view{0zu, N} | std::views::reverse) {
                for(const size_t s : matrix.portrait.shifts_range(i))
                    for(const size_t j_dof : std::ranges::iota_view{0zu, N})
                        if (const size_t j = matrix.portrait.indices[s]; N * j + j_dof > N * i + i_dof)
                            result[i][i_dof] -= matrix.values[s][i_dof][j_dof] * result[j][j_dof];
                result[i][i_dof] = matrix(i, i)[i_dof][i_dof] * result[i][i_dof];
            }
        }

        return result;
    }

public:
    explicit ilu0_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _matrix{std::move(matrix)} {
        if (_matrix.rows() != _matrix.cols())
            throw std::invalid_argument{"ILU0 preconditioner requires a square matrix."};
        validate_sparse_matrix(_matrix);
        compute(_matrix);
    }

    // Solves (L U) x = rhs via forward and backward substitution.
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        if (rhs.size() != matrix().cols())
            throw std::invalid_argument{"ILU0 preconditioner requires rhs vector of the same size as the matrix."};
        return solve(_matrix, rhs);
    }

    const metamath::linear::sparse_matrix<T, I, J>& matrix() const noexcept {
        return _matrix;
    };
};

}

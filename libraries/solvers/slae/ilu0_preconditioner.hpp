#pragma once

#include "preconditioner_base.hpp"

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/types/traits.hpp>
#include <metamath/utils/operators.hpp>

#include <vector>

namespace nonlocal::slae {

// ILU0 (Incomplete LU with zero fill-in) preconditioner for general sparse matrices.
// T may be scalar or square_matrix<floating_point, N> for block structure.
// Factorization: A ≈ L U, sparsity pattern preserved.
// L is stored in the lower triangle (unit diagonal implicit), U in the upper triangle including the diagonal.
template<class T, std::integral I, std::integral J>
class ilu0_preconditioner final : public preconditioner_base<T> {
    metamath::linear::sparse_matrix<T, I, J> _matrix;

    // Computes ILU0 factorization in-place on the sparsity pattern of the matrix.
    void compute() {
        for(const size_t i : std::ranges::iota_view{0zu, _matrix.cols()}) {
            // Process lower entries of row i: for each k < i with (i,k) in pattern
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t k = _matrix.portrait.indices[s]; k < i) {
                    using namespace metamath::linear;
                    // l[i,k] = a[i,k] * inv(u[k,k])
                    _matrix.values[s] *= _matrix(k, k);
                    const T& lik = _matrix.values[s];
                    // Update all remaining entries in row i after position k
                    for(const size_t si : _matrix.portrait.shifts_range(i))
                        if (const size_t j = _matrix.portrait.indices[si]; j > k && _matrix.portrait.contains(k, j))
                            _matrix.values[si] -= lik * _matrix(k, j);
                }
            _matrix(i, i) = metamath::linear::inverse(_matrix(i, i));
        }
    }

public:
    using typename preconditioner_base<T>::entity_t;

    explicit ilu0_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _matrix{std::move(matrix)} {
        if (_matrix.rows() != _matrix.cols())
            throw std::invalid_argument{"ILU0 preconditioner requires a square matrix."};
        compute();
    }

    // Solves (L U) x = rhs via forward and backward substitution.
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        using metamath::operators::operator-=;

        if (rhs.size() != _matrix.cols())
            throw std::invalid_argument{"ILU0 preconditioner requires rhs vector of the same size as the matrix."};
        std::vector<entity_t> result = rhs;

        // Forward substitution: L z = rhs  (L has implicit unit diagonal)
        // z[i] = rhs[i] - sum_{k < i, (i,k) in pattern} l[i,k] * z[k]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t k = _matrix.portrait.indices[s]; k < i)
                    result[i] -= _matrix.values[s] * result[k];

        // Backward substitution: U x = z
        // x[i] = inv(u[i,i]) * (z[i] - sum_{j > i, (i,j) in pattern} u[i,j] * x[j])
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()} | std::views::reverse) {
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i)
                    result[i] -= _matrix.values[s] * result[j];
            result[i] = _matrix(i, i) * result[i];
        }

        return result;
    }
};

}

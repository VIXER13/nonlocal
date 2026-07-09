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
class ilu0_preconditioner final : public preconditioner_base<T, I, J> {
    std::vector<T> _values;
    std::vector<T> _inv_diag;
    const metamath::linear::sparse_matrix_portrait<I, J>* _portrait = nullptr;

    static T invert(const T& val) {
        if constexpr (metamath::types::is_array_v<T>)
            return metamath::linear::inverse(val);
        else
            return T{1} / val;
    }

public:
    using typename preconditioner_base<T, I, J>::entity_t;

    // Computes ILU0 factorization in-place on the sparsity pattern of the matrix.
    void compute(const metamath::linear::sparse_matrix<T, I, J>& matrix) override {
        using namespace metamath::linear;
        const size_t n = matrix.rows();
        _portrait = &matrix.portrait;
        _values = matrix.values;
        _inv_diag.resize(n);

        for (size_t i = 0; i < n; ++i) {
            // Process lower entries of row i: for each k < i with (i,k) in pattern
            for (J s = matrix.portrait.shifts[i]; s < matrix.portrait.shifts[i + 1]; ++s) {
                const size_t k = matrix.portrait.indices[s];
                if (k >= i)
                    break; // entries are sorted; once we reach diagonal, stop

                // l[i,k] = a[i,k] * inv(u[k,k])
                _values[s] *= _inv_diag[k];
                const T& lik = _values[s];

                // Update all remaining entries in row i after position k
                for (J si = matrix.portrait.shifts[i]; si < matrix.portrait.shifts[i + 1]; ++si) {
                    const size_t j = matrix.portrait.indices[si];
                    if (j <= k)
                        continue;
                    if (matrix.portrait.contains(k, j))
                        _values[si] -= lik * _values[matrix.portrait.shift(k, j)];
                }
            }

            _inv_diag[i] = invert(_values[matrix.portrait.shift(i, i)]);
        }
    }

    // Solves (L U) x = rhs via forward and backward substitution.
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        using metamath::operators::operator-=;
        const size_t n = rhs.size();
        std::vector<entity_t> result = rhs;

        // Forward substitution: L z = rhs  (L has implicit unit diagonal)
        // z[i] = rhs[i] - sum_{k < i, (i,k) in pattern} l[i,k] * z[k]
        for (size_t i = 0; i < n; ++i)
            for (J s = _portrait->shifts[i]; s < _portrait->shifts[i + 1]; ++s) {
                const size_t k = _portrait->indices[s];
                if (k >= i)
                    break;
                result[i] -= _values[s] * result[k];
            }

        // Backward substitution: U x = z
        // x[i] = inv(u[i,i]) * (z[i] - sum_{j > i, (i,j) in pattern} u[i,j] * x[j])
        for (size_t i = n; i-- > 0;) {
            for (J s = _portrait->shifts[i]; s < _portrait->shifts[i + 1]; ++s) {
                const size_t j = _portrait->indices[s];
                if (j > i)
                    result[i] -= _values[s] * result[j];
            }
            result[i] = _inv_diag[i] * result[i];
        }

        return result;
    }
};

}

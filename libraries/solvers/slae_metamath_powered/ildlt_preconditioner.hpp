#pragma once

#include "preconditioner_base.hpp"

#include <metamath/linear/fixed_matrix.hpp>
#include <metamath/types/traits.hpp>
#include <metamath/utils/operators.hpp>

#include <vector>

namespace nonlocal::slae {

// Incomplete LDLT factorization preconditioner for symmetric sparse matrices stored in upper triangular format.
// T may be scalar or square_matrix<floating_point, N> for block structure.
// Factorization: A ≈ L D L^T (ILU0 - zero fill-in, sparsity pattern preserved).
// Upper entries store U_{ij} = L_{ji}^T (transpose of lower factor L).
// Diagonal stores D_i (block diagonal factors).
template<class T, std::integral I, std::integral J>
class ildlt_preconditioner final : public preconditioner_base<T, I, J> {
    std::vector<T> _values;
    std::vector<T> _d_values;
    std::vector<T> _inv_d;
    std::vector<std::vector<size_t>> _col_preds;
    const metamath::linear::sparse_matrix_portrait<I, J>* _portrait = nullptr;

public:
    using typename preconditioner_base<T, I, J>::entity_t;

    // Computes incomplete LDLT factorization of the upper symmetric sparse matrix.
    // Sparsity pattern is preserved (ILU0: zero fill-in).
    void compute(const metamath::linear::sparse_matrix<T, I, J>& matrix) override {
        using namespace metamath::linear;
        const size_t n = matrix.rows();
        _portrait = &matrix.portrait;
        _values = matrix.values;
        _d_values.resize(n);
        _inv_d.resize(n);

        // col_preds[j] = rows k < j such that (k, j) is in the upper pattern
        _col_preds.assign(n, {});
        for (size_t k = 0; k < n; ++k)
            for (J s = matrix.portrait.shifts[k]; s < matrix.portrait.shifts[k + 1]; ++s)
                if (const size_t j = matrix.portrait.indices[s]; j > k)
                    _col_preds[j].push_back(k);

        for (size_t i = 0; i < n; ++i) {
            // D_i = A_{ii} - sum_{k in col_preds[i]} transpose(U_{ki}) * D_k * U_{ki}
            T diag = _values[matrix.portrait.shift(i, i)];
            for (const size_t k : _col_preds[i]) {
                const T& uki = _values[matrix.portrait.shift(k, i)];
                diag -= transpose(uki) * _d_values[k] * uki;
            }
            _d_values[i] = diag;
            _inv_d[i] = metamath::linear::inverse(diag);

            // U_{ij} = inv(D_i) * (A_{ij} - sum_{k in col_preds[i] with (k,j) in pattern} transpose(U_{ki}) * D_k * U_{kj})
            for (J s = matrix.portrait.shifts[i]; s < matrix.portrait.shifts[i + 1]; ++s) {
                const size_t j = matrix.portrait.indices[s];
                if (j <= i)
                    continue;
                T val = _values[s];
                for (const size_t k : _col_preds[i])
                    if (matrix.portrait.contains(k, j)) {
                        const T& uki = _values[matrix.portrait.shift(k, i)];
                        const T& ukj = _values[matrix.portrait.shift(k, j)];
                        val -= transpose(uki) * _d_values[k] * ukj;
                    }
                _values[s] = _inv_d[i] * val;
            }
        }
    }

    // Solves (L D L^T) x = rhs via three triangular sweeps.
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        using metamath::operators::operator-=;
        const size_t n = rhs.size();
        std::vector<entity_t> result = rhs;

        // Forward substitution: L y = rhs
        // y[i] = rhs[i] - sum_{k in col_preds[i]} L_{ik} * y[k]
        //      = rhs[i] - sum_{k in col_preds[i]} transpose(U_{ki}) * y[k]
        for (size_t i = 0; i < n; ++i)
            for (const size_t k : _col_preds[i])
                result[i] -= transpose(_values[_portrait->shift(k, i)]) * result[k];

        // Diagonal solve: D z = y  ->  z[i] = inv(D_i) * y[i]
        for (size_t i = 0; i < n; ++i)
            result[i] = _inv_d[i] * result[i];

        // Backward substitution: L^T x = z
        // x[i] = z[i] - sum_{j > i, (i,j) in upper pattern} U_{ij} * x[j]
        for (size_t i = n; i-- > 0;)
            for (J s = _portrait->shifts[i]; s < _portrait->shifts[i + 1]; ++s) {
                const size_t j = _portrait->indices[s];
                if (j > i)
                    result[i] -= _values[s] * result[j];
            }

        return result;
    }
};

}

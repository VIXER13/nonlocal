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
class ildlt_preconditioner final : public preconditioner_base<T> {
    metamath::linear::sparse_matrix<T, I, J> _matrix;
    std::vector<T> _inv_d;
    std::vector<std::vector<size_t>> _col_preds;

    // Computes incomplete LDLT factorization of the upper symmetric sparse matrix.
    // Sparsity pattern is preserved (ILU0: zero fill-in).
    void compute() {
        using namespace metamath::linear;
        const size_t n = _matrix.cols();
        // col_preds[j] = rows k < j such that (k, j) is in the upper pattern
        _col_preds.assign(n, {});
        for(const size_t k : std::ranges::iota_view{0zu, n})
            for(const size_t s : _matrix.portrait.shifts_range(k))
                if (const size_t j = _matrix.portrait.indices[s]; j > k)
                    _col_preds[j].push_back(k);

        std::vector<T> d_values(n);
        for(const size_t i : std::ranges::iota_view{0zu, n}) {
            // D_i = A_{ii} - sum_{k in col_preds[i]} transpose(U_{ki}) * D_k * U_{ki}
            d_values[i] = _matrix(i, i);
            for (const size_t k : _col_preds[i]) {
                const T& uki = _matrix(k, i);
                d_values[i] -= transpose(uki) * d_values[k] * uki;
            }
            _inv_d[i] = metamath::linear::inverse(d_values[i]);

            // U_{ij} = inv(D_i) * (A_{ij} - sum_{k in col_preds[i] with (k,j) in pattern} transpose(U_{ki}) * D_k * U_{kj})
            for (const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i) {
                    for (const size_t k : _col_preds[i])
                        if (_matrix.portrait.contains(k, j))
                            _matrix.values[s] -= transpose(_matrix(k, i)) * d_values[k] * _matrix(k, j);
                    _matrix.values[s] = _inv_d[i] * _matrix.values[s];
                }
        }
    }

public:
    using typename preconditioner_base<T>::entity_t;

    explicit ildlt_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _matrix{std::move(matrix)}, _inv_d(_matrix.cols()) {
        if (_matrix.rows() != _matrix.cols())
            throw std::invalid_argument{"Incomplete LDLT preconditioner requires a square matrix."};
        compute();
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
        for(const size_t i : std::ranges::iota_view{0zu, n})
            for (const size_t k : _col_preds[i])
                result[i] -= transpose(_matrix(k, i)) * result[k];

        // Diagonal solve: D z = y  ->  z[i] = inv(D_i) * y[i]
        for (const size_t i : std::ranges::iota_view{0zu, n})
            result[i] = _inv_d[i] * result[i];

        // Backward substitution: L^T x = z
        // x[i] = z[i] - sum_{j > i, (i,j) in upper pattern} U_{ij} * x[j]
        for (const size_t i : std::ranges::iota_view{0zu, n} | std::views::reverse)
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i)
                    result[i] -= _matrix.values[s] * result[j];

        return result;
    }
};

}

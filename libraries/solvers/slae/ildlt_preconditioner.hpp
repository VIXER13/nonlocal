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
//
// Storage format for diagonal blocks after factorization:
//   scalar T:  stored as inv(D_i)
//   block T:   stored[r][r] = inv(d_r),  stored[r][c>r] = L[c][r],  stored[r][c<r] = L[r][c]
//              where d_r are scalar diagonal pivots from internal LDLT of D_i.
// Off-diagonal upper blocks store U_{ij} factors.

namespace {

// Factorizes the symmetric diagonal block in-place as LDLT.
// Scalar: stores inv(d).
// Block: row-by-row LDLT within the block.
template<std::floating_point T>
void factorize_diagonal_ldlt(T& stored, const T& full_diag) noexcept {
    stored = T{1} / full_diag;
}

template<std::floating_point T, size_t N>
void factorize_diagonal_ldlt(
    metamath::linear::square_matrix<T, N>& stored,
    const metamath::linear::square_matrix<T, N>& full_diag) noexcept {
    stored = full_diag;
    for (size_t r = 0; r < N; ++r) {
        const T dr = stored[r][r];
        stored[r][r] = T{1} / dr;  // store inv(d_r)
        for (size_t r2 = r + 1; r2 < N; ++r2) {
            stored[r2][r] /= dr;           // L[r2][r]
            stored[r][r2] = stored[r2][r]; // L^T: symmetric storage
        }
        // Schur complement (symmetric)
        for (size_t r2 = r + 1; r2 < N; ++r2)
            for (size_t c = r + 1; c < N; ++c)
                stored[r2][c] -= stored[r2][r] * dr * stored[c][r];
    }
}

// Computes inv(D_i) * B column by column using the stored LDLT factored diagonal.
// Scalar: direct multiply.
// Block: forward sub (L), diagonal scale (inv(D_diag)), backward sub (L^T).
template<std::floating_point T>
T apply_inv_diag_left(const T& stored_diag, const T& block) noexcept {
    return stored_diag * block;
}

template<std::floating_point T, size_t N>
metamath::linear::square_matrix<T, N> apply_inv_diag_left(
    const metamath::linear::square_matrix<T, N>& stored_diag,
    const metamath::linear::square_matrix<T, N>& B) noexcept {
    // Solve D_i * X = B where D_i = L * D_diag * L^T (factored in stored_diag)
    // Process column by column.
    metamath::linear::square_matrix<T, N> X = B;
    for (size_t col = 0; col < N; ++col) {
        // Forward sub: inv(L) * B[:,col]
        for (size_t r = 1; r < N; ++r)
            for (size_t c = 0; c < r; ++c)
                X[r][col] -= stored_diag[r][c] * X[c][col];
        // Diagonal scale: inv(D_diag)
        for (size_t r = 0; r < N; ++r)
            X[r][col] *= stored_diag[r][r];
        // Backward sub: inv(L^T)
        for (size_t r = N; r-- > 0;)
            for (size_t c = r + 1; c < N; ++c)
                X[r][col] -= stored_diag[r][c] * X[c][col];
    }
    return X;
}

// Applies inv(D_i) to a vector using the stored LDLT factored diagonal.
// Scalar: x *= inv(d).
// Block: forward sub (L), diagonal scale, backward sub (L^T).
template<std::floating_point T>
void apply_inv_diag_vec(T& x, const T& stored_diag) noexcept {
    x *= stored_diag;
}

template<std::floating_point T, size_t N>
void apply_inv_diag_vec(
    std::array<T, N>& x,
    const metamath::linear::square_matrix<T, N>& stored_diag) noexcept {
    // Forward sub: inv(L)
    for (size_t r = 1; r < N; ++r)
        for (size_t c = 0; c < r; ++c)
            x[r] -= stored_diag[r][c] * x[c];
    // Diagonal scale: inv(D_diag)
    for (size_t r = 0; r < N; ++r)
        x[r] *= stored_diag[r][r];
    // Backward sub: inv(L^T)
    for (size_t r = N; r-- > 0;)
        for (size_t c = r + 1; c < N; ++c)
            x[r] -= stored_diag[r][c] * x[c];
}

} // anonymous namespace

template<class T, std::integral I, std::integral J>
class ildlt_preconditioner final : public preconditioner_base<T> {
    metamath::linear::sparse_matrix<T, I, J> _matrix;
    std::vector<std::vector<size_t>> _col_preds; // TODO: consider using metamath::linear::sparse_matrix_portrait for column predecessors
                                                 // or reuse the existing portrait of the matrix to avoid extra memory allocation.

    // Computes incomplete LDLT factorization of the upper symmetric sparse matrix.
    // Sparsity pattern is preserved (ILU0: zero fill-in).
    void compute() {
        using namespace metamath::linear;
        // col_preds[j] = rows k < j such that (k, j) is in the upper pattern
        _col_preds.assign(_matrix.cols(), {});
        for(const size_t k : std::ranges::iota_view{0zu, _matrix.cols()})
            for(const size_t s : _matrix.portrait.shifts_range(k))
                if (const size_t j = _matrix.portrait.indices[s]; j > k)
                    _col_preds[j].push_back(k);

        std::vector<T> diagonal(_matrix.cols());
        for(const size_t i : std::ranges::iota_view{0zu, _matrix.cols()}) {
            // D_i = A_{ii} - sum_{k in col_preds[i]} transpose(U_{ki}) * D_k * U_{ki}
            auto& dii = _matrix(i, i);
            diagonal[i] = dii;
            for(const size_t k : _col_preds[i]) {
                const T& uki = _matrix(k, i);
                diagonal[i] -= transpose(uki) * diagonal[k] * uki;
            }
            // Factorize D_i internally as LDLT and store in diagonal position.
            factorize_diagonal_ldlt(dii, diagonal[i]);

            // U_{ij} = inv(D_i) * (A_{ij} - sum_{k in col_preds[i] with (k,j) in pattern} transpose(U_{ki}) * D_k * U_{kj})
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i) {
                    for(const size_t k : _col_preds[i])
                        if (_matrix.portrait.contains(k, j))
                            _matrix.values[s] -= transpose(_matrix(k, i)) * diagonal[k] * _matrix(k, j);
                    _matrix.values[s] = apply_inv_diag_left(dii, _matrix.values[s]);
                }
        }
    }

public:
    using typename preconditioner_base<T>::entity_t;

    explicit ildlt_preconditioner(metamath::linear::sparse_matrix<T, I, J>&& matrix)
        : _matrix{std::move(matrix)} {
        if (_matrix.rows() != _matrix.cols())
            throw std::invalid_argument{"Incomplete LDLT preconditioner requires a square matrix."};
        compute();
    }

    // Solves (L D L^T) x = rhs via three triangular sweeps.
    std::vector<entity_t> solve(const std::vector<entity_t>& rhs) const override {
        using namespace metamath::linear;
        using metamath::operators::operator-=;
        if (rhs.size() != _matrix.cols())
            throw std::invalid_argument{"Incomplete LDLT preconditioner requires rhs vector of the same size as the matrix."};
        std::vector<entity_t> result = rhs;

        // Forward substitution: L y = rhs
        // y[i] = rhs[i] - sum_{k in col_preds[i]} L_{ik} * y[k]
        //      = rhs[i] - sum_{k in col_preds[i]} transpose(U_{ki}) * y[k]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()})
           for (const size_t k : _col_preds[i])
               result[i] -= transpose(_matrix(k, i)) * result[k];

        // Diagonal solve: D z = y  ->  z[i] = inv(D_i) * y[i]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()})
            apply_inv_diag_vec(result[i], _matrix(i, i));

        // Backward substitution: L^T x = z
        // x[i] = z[i] - sum_{j > i, (i,j) in upper pattern} U_{ij} * x[j]
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()} | std::views::reverse)
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i)
                    result[i] -= _matrix.values[s] * result[j];

        return result;
    }

    metamath::linear::sparse_matrix<T, I, J> matrix() const noexcept {
        return _matrix;
    }
};

}

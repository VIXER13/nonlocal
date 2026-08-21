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
//
// Storage format for diagonal blocks after factorization:
//   scalar T:  stored as inv(U[i,i])
//   block T:   stored[r][r] = inv(U[r][r]),  stored[r][c>r] = U[r][c],  stored[r][c<r] = L[r][c]
// Off-diagonal lower blocks store the L factors; upper blocks store the U factors.

namespace {

// Computes A * U^{-1} where the stored diagonal encodes the factored format above.
// Scalar: direct multiply (stored = 1/u).
// Block: solves X*U = A column by column left-to-right.
template<std::floating_point T>
T apply_right_inv_upper(const T& A, const T& stored_diag) noexcept {
    return A * stored_diag;
}

template<std::floating_point T, size_t N>
metamath::linear::square_matrix<T, N> apply_right_inv_upper(
    const metamath::linear::square_matrix<T, N>& A,
    const metamath::linear::square_matrix<T, N>& stored_diag) noexcept {
    // Solves X*U = A column by column left-to-right.
    // For column c: X[r][c] = (A[r][c] - sum_{k<c} X[r][k]*U[k][c]) * inv(U[c][c])
    // U[k][c] = stored_diag[k][c] (k < c, upper triangle).
    metamath::linear::square_matrix<T, N> X{};
    for (size_t c = 0; c < N; ++c)
        for (size_t r = 0; r < N; ++r) {
            X[r][c] = A[r][c];
            for (size_t k = 0; k < c; ++k)
                X[r][c] -= X[r][k] * stored_diag[k][c];  // U[k][c] from upper triangle
            X[r][c] *= stored_diag[c][c];                 // multiply by inv(U[c][c])
        }
    return X;
}

// Applies intra-block L forward substitution: x[r] -= sum_{c<r} L[r][c] * x[c].
// Scalar: no-op (no off-diagonal within a scalar "block").
template<std::floating_point T>
void apply_intra_fwd_sub(T&, const T&) noexcept {}

template<std::floating_point T, size_t N>
void apply_intra_fwd_sub(
    std::array<T, N>& x,
    const metamath::linear::square_matrix<T, N>& stored_diag) noexcept {
    for (size_t r = 1; r < N; ++r)
        for (size_t c = 0; c < r; ++c)
            x[r] -= stored_diag[r][c] * x[c];
}

// Applies intra-block U backward substitution: x[r] = inv(U[r][r]) * (x[r] - sum_{c>r} U[r][c]*x[c]).
// Scalar: x *= inv(u) = stored_diag * x.
template<std::floating_point T>
void apply_intra_bwd_sub(T& x, const T& stored_diag) noexcept {
    x *= stored_diag;
}

template<std::floating_point T, size_t N>
void apply_intra_bwd_sub(
    std::array<T, N>& x,
    const metamath::linear::square_matrix<T, N>& stored_diag) noexcept {
    for (size_t r = N; r-- > 0;) {
        for (size_t c = r + 1; c < N; ++c)
            x[r] -= stored_diag[r][c] * x[c];  // U[r][c] from upper triangle
        x[r] *= stored_diag[r][r];              // scale by inv(U[r][r])
    }
}

} // anonymous namespace

template<class T, std::integral I, std::integral J>
class ilu0_preconditioner final : public preconditioner_base<T> {
    metamath::linear::sparse_matrix<T, I, J> _matrix;

    // Factorizes the diagonal block at block row i row by row.
    // For scalar T: inverts the single diagonal element (same as before).
    // For block T: performs row-by-row LU inside the 2D block, storing L factors in the lower
    // triangle and inv(U[r][r]) on the diagonal.  Also propagates the inner L factors to
    // every upper off-diagonal block in row i so that the U blocks are fully updated.
    void factorize_diagonal(const size_t i) {
        if constexpr (metamath::types::is_array_v<T>) {
            using namespace metamath::linear;
            constexpr size_t N = std::tuple_size_v<T>;
            auto& diag = _matrix(i, i);
            for (size_t r = 0; r < N; ++r) {
                diag[r][r] = inverse(diag[r][r]);   // store inv(U[r][r])
                for (size_t r2 = r + 1; r2 < N; ++r2) {
                    const auto L_r2r = diag[r2][r] * diag[r][r];  // L[r2][r] = A[r2][r]/U[r][r]
                    diag[r2][r] = L_r2r;                           // store L factor in lower triangle
                    for (size_t c = r + 1; c < N; ++c)
                        diag[r2][c] -= L_r2r * diag[r][c];        // Schur complement within block
                    // Propagate inner L factor to every upper off-diagonal block in this block row.
                    for (const size_t s : _matrix.portrait.shifts_range(i))
                        if (_matrix.portrait.indices[s] > i)
                            for (size_t c = 0; c < N; ++c)
                                _matrix.values[s][r2][c] -= L_r2r * _matrix.values[s][r][c];
                }
            }
        } else {
            _matrix(i, i) = metamath::linear::inverse(_matrix(i, i));
        }
    }

    // Computes ILU0 factorization in-place on the sparsity pattern of the matrix.
    void compute() {
        using namespace metamath::linear;
        for(const size_t i : std::ranges::iota_view{0zu, _matrix.cols()}) {
            // Process lower entries of row i: for each k < i with (i,k) in pattern.
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t k = _matrix.portrait.indices[s]; k < i) {
                    // l[i,k] = a[i,k] * U[k,k]^{-1}  (right-multiply by upper triangular inverse).
                    // For scalar: same as the old *= stored_inv.
                    // For block: backward column substitution to avoid mixing L and U^{-1} factors.
                    _matrix.values[s] = apply_right_inv_upper(_matrix.values[s], _matrix(k, k));
                    const T& lik = _matrix.values[s];
                    // Update all remaining entries in row i after column k.
                    for(const size_t si : _matrix.portrait.shifts_range(i))
                        if (const size_t j = _matrix.portrait.indices[si]; j > k && _matrix.portrait.contains(k, j))
                            _matrix.values[si] -= lik * _matrix(k, j);
                }
            // Factorize the diagonal block row by row, then update upper off-diagonal blocks.
            factorize_diagonal(i);
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

        // Forward substitution: L z = rhs  (L has implicit unit block-diagonal)
        // z[i] = rhs[i] - sum_{k < i, (i,k) in pattern} L[i,k] * z[k]
        //      followed by intra-block L forward sub for block T.
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()}) {
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t k = _matrix.portrait.indices[s]; k < i)
                    result[i] -= _matrix.values[s] * result[k];
            // For block T: apply lower-triangular L factors stored inside the diagonal block.
            // For scalar T: no-op.
            apply_intra_fwd_sub(result[i], _matrix(i, i));
        }

        // Backward substitution: U x = z
        // x[i] = inv(U[i,i]) * (z[i] - sum_{j > i, (i,j) in pattern} U[i,j] * x[j])
        for(const size_t i : std::ranges::iota_view{0zu, rhs.size()} | std::views::reverse) {
            for(const size_t s : _matrix.portrait.shifts_range(i))
                if (const size_t j = _matrix.portrait.indices[s]; j > i)
                    result[i] -= _matrix.values[s] * result[j];
            // For scalar T: result[i] *= inv(u_ii).
            // For block T: backward sub using stored U entries and inv(U[r][r]) values.
            apply_intra_bwd_sub(result[i], _matrix(i, i));
        }

        return result;
    }
};

}

#pragma once

#include "sparse_matrix.hpp"

#include <metamath/utils/operators.hpp>

namespace metamath::linear {

template<class T, class U, std::integral I, std::integral J>
sparse_matrix<T, I, J>& operator*=(sparse_matrix<T, I, J>& matrix, const U scalar) {
    for(T& value : matrix.values)
         value *= scalar;
     return matrix;
}

template<class T, class U, std::integral I, std::integral J>
sparse_matrix<T, I, J>& operator/=(sparse_matrix<T, I, J>& matrix, const U scalar) {
    for(T& value : matrix.values)
        value /= scalar;
    return matrix;
}

template<class T, class U, std::integral I, std::integral J>
std::vector<U> operator*(const sparse_matrix<T, I, J>& matrix, const std::vector<U>& vector) {
    if (matrix.cols() != vector.size())
        throw std::invalid_argument{"Matrix columns count must match vector size for multiplication."};

    using metamath::operators::operator+=;
    std::vector<U> result(matrix.rows(), U{});
#pragma omp parallel for
    for(size_t row = 0; row < matrix.rows(); ++row)
        for(const size_t shift : matrix.portrait.shifts_range(row))
            result[row] += matrix.values[shift] * vector[matrix.portrait.indices[shift]];
    return result;
}

template<matrix_part Part, class T, class U, std::integral I, std::integral J>
std::vector<U> operator*(const self_adjoint_view<Part, T, I, J>& view, const std::vector<U>& vector) {
    if (view.matrix.cols() != vector.size())
        throw std::invalid_argument{"Matrix columns count must match vector size for multiplication."};

    using metamath::operators::operator+=;
    static constexpr std::conditional_t<Part == matrix_part::Upper, std::greater<>, std::less<>> comparator{};
    std::vector<U> result(view.matrix.rows(), U{});
    for(const size_t row : std::ranges::iota_view{0zu, view.matrix.rows()})
        for(const size_t shift : view.matrix.portrait.shifts_range(row))
            if (const size_t col = view.matrix.portrait.indices[shift]; row == col)
                result[row] += self_adjoint<Part>(view.matrix.values[shift]) * vector[col];
            else if (comparator(col, row)) {
                result[row] += view.matrix.values[shift] * vector[col];
                result[col] += transpose(view.matrix.values[shift]) * vector[row];
            }
    return result;
}

template<class T, std::integral I, std::integral J>
sparse_matrix<T, I, J>& operator+=(sparse_matrix<T, I, J>& lhs, const sparse_matrix<T, I, J>& rhs) {
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
        throw std::invalid_argument{"Matrices must have the same dimensions for addition."};

    for(size_t row = 0; row < lhs.rows(); ++row) {
        size_t lhs_shift = lhs.portrait.shifts[row];
        size_t rhs_shift = rhs.portrait.shifts[row];
        while (lhs_shift < lhs.portrait.shifts[row + 1] && rhs_shift < rhs.portrait.shifts[row + 1]) {
            const size_t lhs_col = lhs.portrait.indices[lhs_shift];
            const size_t rhs_col = rhs.portrait.indices[rhs_shift];
            if (lhs_col == rhs_col) {
                // both matrices have a non-zero element at the same row and col
                lhs.values[lhs_shift] += rhs.values[rhs_shift];
                ++lhs_shift;
                ++rhs_shift;
            } else if (lhs_col < rhs_col) {
                // Move to the next non-zero element in lhs, in rhs it absent.
                ++lhs_shift;
            } else {
                // Insert non-zero element from rhs into lhs
                lhs.portrait.indices.insert(lhs.portrait.indices.begin() + lhs_shift, rhs_col);
                lhs.values.insert(lhs.values.begin() + lhs_shift, rhs.values[rhs_shift]);
                ++lhs_shift;
                ++rhs_shift;
                for(size_t i = row + 1; i <= lhs.rows(); ++i)
                    ++lhs.portrait.shifts[i];
            }
        }

        // Append remaining non-zero elements from rhs if any
        while (rhs_shift < rhs.portrait.shifts[row + 1]) {
            const size_t rhs_col = rhs.portrait.indices[rhs_shift];
            lhs.portrait.indices.insert(lhs.portrait.indices.begin() + lhs_shift, rhs_col);
            lhs.values.insert(lhs.values.begin() + lhs_shift, rhs.values[rhs_shift]);
            ++lhs_shift;
            ++rhs_shift;
            for(size_t i = row + 1; i <= lhs.rows(); ++i)
                ++lhs.portrait.shifts[i];
        }
    }
    return lhs;
}

}
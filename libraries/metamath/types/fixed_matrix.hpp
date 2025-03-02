#pragma once

#include <array>
#include <ranges>
#include <cstddef>

namespace metamath::types {

template<class T, size_t Rows, size_t Cols>
using fixed_matrix = std::array<std::array<T, Cols>, Rows>;

template<class T, size_t N>
using square_matrix = fixed_matrix<T, N, N>;

template<class T, size_t Rows, size_t Cols>
fixed_matrix<T, Rows, Cols> make_fixed_matrix(const std::array<T, Rows * Cols>& arr) {
    fixed_matrix<T, Rows, Cols> matrix;
    for(const size_t row : std::ranges::iota_view{0u, Rows})
        for(const size_t col : std::ranges::iota_view{0u, Cols})
            matrix[row][col] = arr[row * Cols + col];
    return matrix;
}

template<class T, size_t N>
square_matrix<T, N> make_square_matrix(const std::array<T, N * N>& arr) {
    return make_fixed_matrix<T, N, N>(arr);
}

template<class T>
bool is_positive(const square_matrix<T, 2u>& matrix) {
    const T determinant = matrix[0][0] * matrix[0][0] - 2 * matrix[0][0] * matrix[1][1] + 
                          matrix[1][1] * matrix[1][1] + 4 * matrix[0][1] * matrix[1][0];
    if (determinant < T{0})
        return false;
    return matrix[0][0] + matrix[1][1] - std::sqrt(determinant) > T{0};
}

}
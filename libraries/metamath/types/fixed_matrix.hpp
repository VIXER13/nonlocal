#ifndef METAMATH_FIXED_MATRIX_HPP
#define METAMATH_FIXED_MATRIX_HPP

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

}

#endif
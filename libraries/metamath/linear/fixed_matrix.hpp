#pragma once

#include <array>
#include <ranges>
#include <cstddef>
#include <cmath>

namespace metamath::linear {

template<std::floating_point T, size_t Rows, size_t Cols>
using fixed_matrix = std::array<std::array<T, Cols>, Rows>;

template<std::floating_point T, size_t N>
using square_matrix = fixed_matrix<T, N, N>;

template<std::floating_point T>
constexpr T transpose(const T value) noexcept {
    return value;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Cols, Rows> transpose(const fixed_matrix<T, Rows, Cols>& matrix) noexcept {
    fixed_matrix<T, Cols, Rows> transposed{};
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            transposed[col][row] = transpose(matrix[row][col]);
    return transposed;
}

template<std::floating_point T, size_t N>
constexpr T determinant(const square_matrix<T, N>& matrix) noexcept {
    if constexpr (N == 1)
        return matrix[0][0];
    else if constexpr (N == 2)
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    else if constexpr (N == 3)
        return matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
               matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
               matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    else
        static_assert(N < 4, "Determinant calculation is only implemented for 1x1, 2x2, and 3x3 matrices.");
}

template<std::floating_point T, size_t N>
constexpr bool is_positive(const square_matrix<T, N>& matrix) noexcept {
    if constexpr (N == 1)
        return matrix[0][0] > T{0};
    else if constexpr (N == 2)
        return matrix[0][0] > T{0} && determinant(matrix) > T{0};
    else if constexpr (N == 3)
        return matrix[0][0] > T{0} &&
               determinant<T, 2>({matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]}) > T{0} &&
               determinant(matrix) > T{0};
    else
        static_assert(N < 4, "Positive definiteness check is only implemented for 1x1, 2x2, and 3x3 matrices.");
}

template<std::floating_point T, size_t N>
constexpr square_matrix<T, N> inverse(const square_matrix<T, N>& matrix) noexcept {
    if constexpr (N == 1)
        return { T{1} / determinant(matrix) };
    else if constexpr (N == 2) {
        const T det = determinant(matrix);
        return { matrix[1][1] / det, -matrix[0][1] / det,
                -matrix[1][0] / det,  matrix[0][0] / det};
    }
    else if constexpr (N == 3) {
        const T det = determinant(matrix);
        return {(matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / det,
                (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / det,
                (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / det,
                (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / det,
                (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / det,
                (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / det,
                (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / det,
                (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / det,
                (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / det};
    }
    else
        static_assert(N < 4, "Matrix inversion is only implemented for 1x1, 2x2, and 3x3 matrices.");
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator+=(fixed_matrix<T, Rows, Cols>& lhs, const fixed_matrix<T, Rows, Cols>& rhs) noexcept {
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            lhs[row][col] += rhs[row][col];
    return lhs;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator-=(fixed_matrix<T, Rows, Cols>& lhs, const fixed_matrix<T, Rows, Cols>& rhs) noexcept {
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            lhs[row][col] -= rhs[row][col];
    return lhs;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr std::array<T, Rows> operator*(const fixed_matrix<T, Rows, Cols>& matrix, const std::array<T, Cols>& vector) noexcept {
    std::array<T, Rows> result{};
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            result[row] += matrix[row][col] * vector[col];
    return result;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr std::array<T, Cols> operator*(const std::array<T, Rows>& vector, const fixed_matrix<T, Rows, Cols>& matrix) noexcept {
    std::array<T, Cols> result{};
    for (const size_t col : std::ranges::iota_view{0zu, Cols})
        for (const size_t row : std::ranges::iota_view{0zu, Rows})
            result[col] += vector[row] * matrix[row][col];
    return result;
}

template<std::floating_point T, size_t Rows, size_t K, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols> operator*(const fixed_matrix<T, Rows, K>& lhs, const fixed_matrix<T, K, Cols>& rhs) noexcept {
    fixed_matrix<T, Rows, Cols> result{};
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            for (const size_t k : std::ranges::iota_view{0zu, K})
                result[row][col] += lhs[row][k] * rhs[k][col];
    return result;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols>& operator*=(fixed_matrix<T, Rows, Cols>& lhs, const fixed_matrix<T, Cols, Cols>& rhs) noexcept {
    lhs = lhs * rhs;
    return lhs;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols>& operator*=(fixed_matrix<T, Rows, Cols>& lhs, const T value) noexcept {
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            lhs[row][col] *= value;
    return lhs;
}

template<std::floating_point T, size_t N>
constexpr fixed_matrix<T, N, N> operator/(const fixed_matrix<T, N, N>& lhs, const fixed_matrix<T, N, N>& rhs) noexcept {
    return lhs * inverse(rhs);
}

template<std::floating_point T, size_t N>
constexpr fixed_matrix<T, N, N>& operator/=(fixed_matrix<T, N, N>& lhs, const fixed_matrix<T, N, N>& rhs) noexcept {
    lhs *= inverse(rhs);
    return lhs;
}

template<std::floating_point T, size_t Rows, size_t Cols>
constexpr fixed_matrix<T, Rows, Cols>& operator/=(fixed_matrix<T, Rows, Cols>& lhs, const T value) noexcept {
    for (const size_t row : std::ranges::iota_view{0zu, Rows})
        for (const size_t col : std::ranges::iota_view{0zu, Cols})
            lhs[row][col] /= value;
    return lhs;
}

}
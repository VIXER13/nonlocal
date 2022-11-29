#ifndef METAMATH_FIXED_MATRIX_HPP
#define METAMATH_FIXED_MATRIX_HPP

#include <array>
#include <cinttypes>

namespace metamath::types {

template<class T, size_t Rows, size_t Cols>
using fixed_matrix = std::array<std::array<T, Cols>, Rows>;

template<class T, size_t N>
using square_matrix = fixed_matrix<T, N, N>;

}

#endif
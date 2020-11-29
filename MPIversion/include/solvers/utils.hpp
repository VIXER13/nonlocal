#ifndef NONLOCAL_SOLVERS_UTILS_HPP
#define NONLOCAL_SOLVERS_UTILS_HPP

#include <array>
#include "metamath.hpp"

namespace nonlocal::utils {

template<class T, size_t N>
T distance(const std::array<T, N>& A, const std::array<T, N>& B) noexcept {
    T sum = 0;
    for(size_t i = 0; i < N; ++i)
        sum += metamath::function::power<2>(A[i] - B[i]);
    return sqrt(sum);
}

template<class T, size_t N>
std::array<T, N> operator*(std::array<T, N> arr, const T val) noexcept {
    for(size_t i = 0; i < N; ++i)
        arr[i] *= val;
    return arr;
}

template<class T, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& lhs, const std::array<T, N>& rhs) noexcept {
    for(size_t i = 0; i < N; ++i)
        lhs[i] += rhs[i];
    return lhs;
}

}

#endif
#ifndef METAMATH_DISTANCE_HPP
#define METAMATH_DISTANCE_HPP

#include "power.hpp"

#include <array>
#include <cmath>

namespace metamath::functions {

template<class T, size_t N>
constexpr T distance(const std::array<T, N>& lhs, const std::array<T, N>& rhs) noexcept {
    T sum = 0;
    for(size_t i = 0; i < N; ++i)
        sum += metamath::functions::power<2>(lhs[i] - rhs[i]);
    return std::sqrt(sum);
}

}

#endif
#pragma once

#include <cmath>
#include <concepts>

namespace metamath::functions {

#ifdef __clang__
template<std::floating_point T>
constexpr T sqrt(T value) noexcept {
    T result = value;
    for(T last = T{0}; result != last; result = T{0.5} * (result + value / result))
        last = result;
    return result;
}
#else
template<std::floating_point T> 
constexpr T sqrt(T value) noexcept { return std::sqrt(value); }
#endif

}
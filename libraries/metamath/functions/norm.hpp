#pragma once

#include "power.hpp"

#include <metamath/utils/operators.hpp>
#include <metamath/utils/constants.hpp>
#include <metamath/types/traits.hpp>

#include <array>
#include <cmath>
#include <numeric>

namespace metamath::functions {

template<size_t Exp = 2, class T, size_t D>
T powered_norm(const std::array<T, D>& x) {
    return std::accumulate(x.begin(), x.end(), T{0}, [](const T sum, const T x) {
        if constexpr (Exp & 1)
            return sum + power<Exp>(std::abs(x));
        return sum + power<Exp>(x);
    });
}

template<class T, size_t D, types::arithmetic Exp>
T powered_norm(const std::array<T, D>& x, const Exp exp) {
    return std::accumulate(x.begin(), x.end(), T{0}, [exp](const T sum, const T x) {
        return sum + power(std::abs(x), exp);
    });
}

template<size_t Exp = 2, class T, size_t D>
T norm(const std::array<T, D>& x) {
    if constexpr (Exp == 1)
        return powered_norm<Exp>(x);
    if constexpr (Exp == 2)
        return std::sqrt(powered_norm<Exp>(x));
    if constexpr (Exp == 3)
        return std::cbrt(powered_norm<Exp>(x));
    if constexpr (Exp == constants::Infinity<size_t>)
        return std::abs(*std::max_element(x.begin(), x.end(), 
            [](const T a, const T b) { return std::abs(a) < std::abs(b); }));
    return power(powered_norm<Exp>(x), T{1} / Exp);
}

template<class T, size_t D, types::arithmetic Exp>
T norm(const std::array<T, D>& x, const Exp exp) {
    return power(powered_norm(x, exp), T{1} / exp);
}

}
#pragma once

#include <cmath>
#include <concepts>
#include <type_traits>

namespace metamath::functions {

template<auto N, class T, std::enable_if_t<std::is_integral_v<decltype(N)>, bool> = true>
constexpr T power(const T x) noexcept {
    if constexpr (N < 0)
        return 1 / power<-N>(x);
    else if constexpr (N == 0)
        return 1;
    else if constexpr (N == 1)
        return x;
    else if constexpr (N & 1)
        return x * power<N - 1>(x);
    else {
        const T temp = power<N / 2>(x);
        return temp * temp;
    }
}

template<class T, std::integral Exp>
constexpr T power(T x, Exp exp) noexcept {
    if (exp < 0)
        return T{1} / power(x, -exp);
    T result = T{1};
    while (exp > 0) {
        if (exp & 1) {
            result *= x;
            --exp;
        } else {
            x *= x;
            exp >>= 1;
        }
    }
    return result;
}

template<class T, std::floating_point Exp>
constexpr T power(const T x, Exp exp) noexcept {
    return std::pow(x, exp);
}

}
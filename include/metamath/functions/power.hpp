#ifndef METAMATH_POWER_HPP
#define METAMATH_POWER_HPP

#include <cinttypes>

namespace metamath::function {

template<uintmax_t N, class T>
constexpr T power_u(const T& x) noexcept {
    if constexpr (N == 0)
        return 1;
    else if constexpr (N == 1)
        return x;
    else if constexpr (N % 2)
        return x * power_u<N - 1>(x);
    else {
        const T temp = power_u<N / 2>(x);
        return temp * temp;
    }
}

template<uintmax_t N, class T>
constexpr T power_m(const T& x) noexcept {
    return T{1} / power_u<N>(x);
}

template<intmax_t N, class T>
constexpr T power(const T& x) noexcept {
    if constexpr (N >= 0)
        return power_u<N>(x);
    else
        return power_m<uintmax_t{-(N + 1)} + 1>(x);
}

}

#endif
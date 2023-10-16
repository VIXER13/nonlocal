#ifndef METAMATH_POWER_HPP
#define METAMATH_POWER_HPP

#include <type_traits>

namespace metamath::functions {

template<auto N, class T, std::enable_if_t<std::is_integral_v<decltype(N)>, bool> = true>
constexpr T power(const T& x) {
    if constexpr (N < 0)
        return 1 / power<-N>(x);
    else if constexpr (N == 0)
        return 1;
    else if constexpr (N == 1)
        return x;
    else if constexpr (N % 2)
        return x * power<N - 1>(x);
    else {
        const T temp = power<N / 2>(x);
        return temp * temp;
    }
}

}

#endif
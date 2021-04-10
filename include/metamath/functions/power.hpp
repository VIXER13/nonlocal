#ifndef METAMATH_POWER_HPP
#define METAMATH_POWER_HPP

#include <cinttypes>
#include <type_traits>

namespace metamath::function {

template<intmax_t N, class T>
constexpr T power(const T x) {
    if constexpr (N > 0 && N % 2 == 0) {
        const T temp = power<N / 2>(x);
        return temp * temp;
    } else if constexpr (N > 0 && N % 2 == 1)
        return x * power<N - 1>(x);
    else if constexpr (N < 0)
        return 1 / power<-N>(x);
    return 1;
}

};

#endif
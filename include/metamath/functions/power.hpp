#ifndef METAMATH_POWER_HPP
#define METAMATH_POWER_HPP

#include <cinttypes>

namespace metamath::function {

template<intmax_t N, class T>
constexpr T power(const T& x) {
    if constexpr (N == 0)
        return 1;
    else if constexpr (N < 0)
        return 1 / (x * power<-(N + 1)>(x));
    else if constexpr (N % 2)
        return x * power<N - 1>(x);
    else {
        const T temp = power<N / 2>(x);
        return temp * temp;
    }
}

};

#endif
#ifndef METAMATH_POWER_HPP
#define METAMATH_POWER_HPP

#include <cinttypes>
#include <type_traits>

namespace metamath::function {

template<intmax_t N, class Type> constexpr std::enable_if_t<(N > 0 && N % 2 == 0), Type> power(const Type);
template<intmax_t N, class Type> constexpr std::enable_if_t<(N > 0 && N % 2 == 1), Type> power(const Type);

template<intmax_t N, class Type>
constexpr std::enable_if_t<N < 0, Type> power(const Type x) {
    return 1 / power<-N>(x);
}

template<intmax_t N, class Type>
constexpr std::enable_if_t<N == 0, Type> power(const Type) {
    return 1;
}

template<intmax_t N, class Type>
constexpr std::enable_if_t<(N > 0 && N % 2 == 0), Type> power(const Type x) {
    const Type temp = power<N / 2>(x);
    return temp * temp;
}

template<intmax_t N, class Type>
constexpr std::enable_if_t<(N > 0 && N % 2 == 1), Type> power(const Type x) {
    return x * power<N - 1>(x);
}

};

#endif
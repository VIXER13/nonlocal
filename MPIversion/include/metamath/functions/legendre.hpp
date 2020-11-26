#ifndef METAMATH_FUNCTIONS_LEGENDRE_HPP
#define METAMATH_FUNCTIONS_LEGENDRE_HPP

#include <type_traits>
#include <stdexcept>
#include "multinomial.hpp"
#include "power.hpp"

namespace metamath::function {

class _legendre {
    _legendre() = delete;

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<K == 0, T> legendre_term(const T x) {
        return power<N>(x-1);
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K > 0), T> legendre_term(const T x) {
        static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
        return power<2, int64_t>(constants::binomial<N, K>{}) * power<N-K>(x-1) * power<K>(x+1) + legendre_term<N, K-1>(x);
    }

public:
    template<uintmax_t N, class T, bool With_Check>
    friend constexpr T legendre(const T x);
};

template<uintmax_t N, class T, bool With_Check>
constexpr T legendre(const T x) {
    if constexpr (With_Check) {
        if(x < -1 || x > 1)
            throw std::domain_error{"Argument out of range in legendre."};
    }
    return _legendre::legendre_term<N, N>(x) / power<N, T>(2);
}

}

#endif
#ifndef METAMATH_FUNCTIONS_LAGUERRE_HPP
#define METAMATH_FUNCTIONS_LAGUERRE_HPP

#include <type_traits>
#include <stdexcept>
#include "multinomial.hpp"
#include "power.hpp"

namespace metamath::function {

class _laguerre {
    _laguerre() = delete;

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K == 0), T> laguerre_term(const T) {
        return 1;
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K > 0), T> laguerre_term(const T x) {
        static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
        constexpr int64_t coeff = (K % 2 ? -1 : 1) * constants::binomial<N, K>{};
        return coeff * power<K>(x) / constants::factorial<K>{} + laguerre_term<N, K-1>(x);
    }

public:
    template<uintmax_t N, class T, bool With_Check>
    friend constexpr T laguerre(const T x);
};

template<uintmax_t N, class T, bool With_Check>
constexpr T laguerre(const T x) {
    if constexpr (With_Check) {
        if(x < 0)
            throw std::domain_error{"Negative argument in laguerre."};
    }
    return _laguerre::laguerre_term<N, N>(x);
}

}

#endif
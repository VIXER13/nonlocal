#ifndef METAMATH_FUNCTIONS_CHEBYSHEV_HPP
#define METAMATH_FUNCTIONS_CHEBYSHEV_HPP

#include <type_traits>
#include <stdexcept>
#include "multinomial.hpp"
#include "power.hpp"

namespace metamath::function {

class _chebyshev {
    _chebyshev() = delete;

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K > 0), T> chebyshev_first_kind_term(const T x) {
        static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
        return constants::binomial<N, 2*K>{} * power<K>(x*x - 1) * power<N-2*K>(x) + chebyshev_first_kind_term<N, K-1>(x);
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<K == 0, T> chebyshev_first_kind_term(const T x) {
        return power<N>(x);
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K > 0), T> chebyshev_second_kind_term(const T x) {
        static_assert(std::is_floating_point_v<T>, "The T must be floating point.");
        return constants::binomial<N+1, 2*K+1>{} * power<K>(x*x - 1) * power<N-2*K>(x) + chebyshev_second_kind_term<N, K-1>(x);
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<K == 0, T> chebyshev_second_kind_term(const T x) {
        return constants::binomial<N+1, 1>{} * power<N>(x);
    }

public:
    template<uintmax_t N, class T>
    friend constexpr T chebyshev_first_kind(const T x);

    template<uintmax_t N, class T>
    friend constexpr T chebyshev_second_kind(const T x);
};

template<uintmax_t N, class T>
constexpr T chebyshev_first_kind(const T x) {
    if(x < -1 || x > 1)
        throw std::domain_error{"Argument out of range in math_meta::chebyshev_first_kind."};
    return _chebyshev::chebyshev_first_kind_term<N, N/2>(x);
}

template<uintmax_t N, class T>
constexpr T chebyshev_second_kind(const T x) {
    if(x < -1 || x > 1)
        throw std::domain_error{"Argument out of range in math_meta::chebyshev_second_kind."};
    return _chebyshev::chebyshev_second_kind_term<N, N/2>(x);
}

}

#endif
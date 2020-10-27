#ifndef SYMDIFF_INTEGRAL_CONSTANT_HPP
#define SYMDIFF_INTEGRAL_CONSTANT_HPP

#include "expression.hpp"

namespace metamath::symdiff {

template<class T, T N>
struct integral_constant : expression<integral_constant<T, N>>,
                           std::integral_constant<T, N> {
    template<uintmax_t X>
    using derivative_type = integral_constant<int8_t, 0>;

    template<typename U>
    constexpr T operator()(const U&) const {
        return std::integral_constant<T, N>{};
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return derivative_type<X>{};
    }
};

template<class E>
struct is_integral_constant : std::false_type {};

template<class T, T N>
struct is_integral_constant<integral_constant<T, N>> : std::true_type {};

template<class E>
struct integral_constant_v : std::false_type {};

template<class T, T N>
struct integral_constant_v<integral_constant<T, N>> : std::integral_constant<T, N> {};

template<class E>
using integral_constant_t = typename integral_constant_v<E>::value_type;

}

#endif
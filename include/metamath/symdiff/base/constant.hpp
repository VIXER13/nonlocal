#ifndef SYMDIFF_CONSTANT_HPP
#define SYMDIFF_CONSTANT_HPP

#include "integral_constant.hpp"

namespace metamath::symdiff {

template<class T>
struct constant : expression<constant<T>> {
    using value_type = T;
    template<uintmax_t X>
    using derivative_type = integral_constant<intmax_t, 0>;

    const T value = 0;

    constexpr constant(const T& value) : value{value} {}

    template<class U>
    constexpr T operator()(const U&) const {
        return value;
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return derivative_type<X>{};
    }
};

template<class E>
struct is_constant : std::false_type {};

template<class T>
struct is_constant<constant<T>> : std::true_type {};

template<class T, T N>
struct is_constant<integral_constant<T, N>> : std::true_type {};

template<class E>
struct constant_type : std::false_type {};

template<class T>
struct constant_type<constant<T>> : constant<T> {};

template<class T, T N>
struct constant_type<integral_constant<T, N>> : integral_constant<T, N> {};

template<class E>
using constant_t = typename constant_type<E>::value_type;

}

#endif
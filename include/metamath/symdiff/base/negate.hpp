#ifndef SYMDIFF_NEGATE_HPP
#define SYMDIFF_NEGATE_HPP

#include "constant.hpp"

namespace metamath::symdiff {

template<class E>
class negate;

template<class E>
using negate_type = std::conditional_t<
    is_integral_constant<E>{},
    integral_constant<integral_constant_t<E>, -integral_constant_v<E>{}>,
    std::conditional_t<is_constant<E>{}, E, negate<E>>
>;

template<class E>
class negate : public expression<negate<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = negate_type<typename E::template derivative_type<X>>;

    constexpr explicit negate(const expression<E>& e) : e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(-e(x)) {
        return -e(x);
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return -e.template derivative<X>();
    }
};

template<class T, T N>
constexpr integral_constant<T, -N> operator-(const integral_constant<T, N>&) {
    return integral_constant<T, -N>{};
}

template<class T>
constexpr constant<T> operator-(const constant<T>& e) {
    return constant<T>{-e.value};
}

template<class E>
constexpr negate<E> operator-(const expression<E>& e) {
    return negate<E>{e};
}

}

#endif
#ifndef SYMBOLIC_NEGATE_HPP
#define SYMBOLIC_NEGATE_HPP

#include "constant.hpp"
#include "unary_expression.hpp"

namespace metamath::symbolic {

template<class E>
class negate : public unary_expression<E, negate> {
    using _base = unary_expression<E, negate>;

public:
    constexpr explicit negate(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return -_base::expr()(args...);
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return -_base::expr().template derivative<X>();
    }
};

template<class E>
constexpr negate<E> operator-(const expression<E>& e) {
    return negate<E>{e()};
}

template<auto N>
constexpr integral_constant<-N> simplify(const negate<integral_constant<N>>) noexcept {
    return {};
}

template<class T>
constexpr constant<T> simplify(const negate<constant<T>>& c) {
    return constant{-c()};
}

template<class E>
constexpr auto simplify(const negate<negate<E>>& e) {
    return simplify(e.expr().expr());
}

}

#endif
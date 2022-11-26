#ifndef SYMBOLIC_SQRT_EXPRESSION_HPP
#define SYMBOLIC_SQRT_EXPRESSION_HPP

#include "divides.hpp"
#include "integral_constant.hpp"
#include "unary_expression.hpp"

#include <cmath>

namespace SYMBOLIC_NAMESPACE {

template<class E>
class sqrt_expression : public unary_expression<E, sqrt_expression> {
    using _base = unary_expression<E, sqrt_expression>;

public:
    constexpr explicit sqrt_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    auto operator()(const Args&... args) const {
        return std::sqrt(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return _base::expr().template derivative<X>() / (integral_constant<2>{} * sqrt(_base::expr()));
    }
};

template<class E>
constexpr sqrt_expression<E> sqrt(const expression<E>& e) noexcept {
    return sqrt_expression<E>{e()};
}

}

#endif
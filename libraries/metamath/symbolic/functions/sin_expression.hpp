#pragma once

#include "multiplies.hpp"
#include "unary_expression.hpp"

#include <cmath>

namespace metamath::symbolic {

template<class E>
struct cos_expression;

template<class E>
class sin_expression : public unary_expression<E, sin_expression> {
    using _base = unary_expression<E, sin_expression>;

public:
    constexpr explicit sin_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    auto operator()(const Args&... args) const {
        return std::sin(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return _base::expr().template derivative<X>() * cos(_base::expr());
    }
};

template<class E>
constexpr sin_expression<E> sin(const expression<E>& e) noexcept {
    return sin_expression<E>{e()};
}

}
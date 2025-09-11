#pragma once

#include "sign_expression.hpp"

#include <metamath/symbolic/base/multiplies.hpp>

#include <cmath>

namespace metamath::symbolic {

template<class E>
class abs_expression : public unary_expression<E, abs_expression> {
    using _base = unary_expression<E, abs_expression>;

public:
    constexpr explicit abs_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::abs(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return _base::expr().template derivative<X>() * sign(_base::expr());
    }
};

template<class E>
constexpr abs_expression<E> abs(const expression<E>& e) noexcept {
    return abs_expression<E>{e()};
}

}
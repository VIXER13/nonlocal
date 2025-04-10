#pragma once

#include <metamath/symbolic/base/multiplies.hpp>
#include <metamath/symbolic/base/unary_expression.hpp>

#include <cmath>

namespace metamath::symbolic {

template<class E>
class exp_expression : public unary_expression<E, exp_expression> {
    using _base = unary_expression<E, exp_expression>;

public:
    constexpr explicit exp_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::exp(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return _base::expr().template derivative<X>() * exp(_base::expr());
    }
};

template<class E>
constexpr exp_expression<E> exp(const expression<E>& e) noexcept {
    return exp_expression<E>{e()};
}

}
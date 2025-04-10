#pragma once

#include "exp_expression.hpp"

#include <numbers>

namespace metamath::symbolic {

template<class E, auto T>
class erf_expression : public unary_expression<E, erf_expression, T> {
    using _base = unary_expression<E, erf_expression, T>;

public:
    constexpr explicit erf_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::erf(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return integral_constant<2>{} * std::numbers::inv_sqrtpi_v<decltype(T)> *
               _base::expr().template derivative<X>() * exp(-power<2>(_base::expr()));
    }
};

template<class T, class E>
constexpr erf_expression<E, T{}> erf(const expression<E>& e) noexcept {
    return erf_expression<E, T{}>{e()};
}

}
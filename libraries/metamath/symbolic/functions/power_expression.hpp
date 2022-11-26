#ifndef SYMBOLIC_POWER_EXPRESSION_HPP
#define SYMBOLIC_POWER_EXPRESSION_HPP

#include "functions/power.hpp"

#include "integral_constant.hpp"
#include "multiplies.hpp"
#include "unary_expression.hpp"

#include <type_traits>

namespace SYMBOLIC_NAMESPACE {

template<class E, auto N>
class power_expression : public unary_expression<E, power_expression, N> {
    using _base = unary_expression<E, power_expression, N>;

public:
    constexpr explicit power_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return metamath::functions::power<N>(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return integral_constant<N>{} * _base::expr().template derivative<X>() * power<N-1>(_base::expr());
    }
};

template<auto N, class E>
constexpr power_expression<E, N> power(const expression<E>& e) noexcept {
    return power_expression<E, N>{e()};
}

template<class E, auto N, std::enable_if_t<N == 0, bool> = true>
constexpr integral_constant<1> simplify(const power_expression<E, N>&) {
    return {};
}

template<class E, auto N, std::enable_if_t<N == 1, bool> = true>
constexpr auto simplify(const power_expression<E, N>& e) {
    return simplify(e.expr());
}

template<class E, auto N, auto M>
constexpr auto simplify(const power_expression<power_expression<E, M>, N>& e) {
    return power_expression<E, N * M>{simplify(e.expr().expr())};
}

}

#endif
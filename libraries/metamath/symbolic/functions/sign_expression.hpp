#pragma once

#include <metamath/symbolic/base/integral_constant.hpp>
#include <metamath/symbolic/base/unary_expression.hpp>

namespace metamath::symbolic {

template<class E>
class sign_expression : public unary_expression<E, sign_expression> {
    using _base = unary_expression<E, sign_expression>;

public:
    constexpr explicit sign_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        const auto value = _base::expr()(args...);
        return value < 0 ? -1 :
               value > 0 ?  1 : 0;
    }

    template<auto X>
    constexpr integral_constant<0> derivative() const noexcept {
        return {};
    }
};

template<class E>
constexpr sign_expression<E> sign(const expression<E>& e) {
    return sign_expression<E>{e()};
}

}
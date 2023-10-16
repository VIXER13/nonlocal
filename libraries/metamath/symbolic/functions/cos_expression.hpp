#ifndef SYMBOLIC_COS_EXPRESSION_HPP
#define SYMBOLIC_COS_EXPRESSION_HPP

#include "multiplies.hpp"
#include "negate.hpp"

#include <cmath>

namespace metamath::symbolic {

template<class E>
struct sin_expression;

template<class E>
class cos_expression : public unary_expression<E, cos_expression> {
    using _base = unary_expression<E, cos_expression>;

public:
    constexpr explicit cos_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::cos(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return -_base::expr().template derivative<X>() * sin(_base::expr());
    }
};

template<class E>
constexpr cos_expression<E> cos(const expression<E>& e) noexcept {
    return cos_expression<E>{e()};
}

}

#endif
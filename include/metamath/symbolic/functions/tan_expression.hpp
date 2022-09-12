#ifndef SYMBOLIC_TAN_EXPRESSION_HPP
#define SYMBOLIC_TAN_EXPRESSION_HPP

#include "cos_expression.hpp"
#include "power_expression.hpp"

namespace SYMBOLIC_NAMESPACE {

template<class E>
class tan_expression : public unary_expression<E, tan_expression> {
    using _base = unary_expression<E, tan_expression>;

public:
    constexpr explicit tan_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::tan(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const {
        return _base::expr().template derivative<X>() / power<2>(cos(_base::expr()));
    }
};

template<class E>
constexpr tan_expression<E> tan(const expression<E>& e) noexcept {
    return tan_expression<E>{e};
}

}

#endif
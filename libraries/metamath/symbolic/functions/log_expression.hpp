#ifndef SYMBOLIC_LOG_EXPRESSION_HPP
#define SYMBOLIC_LOG_EXPRESSION_HPP

#include "divides.hpp"
#include "unary_expression.hpp"

#include <cmath>

namespace metamath::symbolic {

template<class E>
class log_expression : public unary_expression<E, log_expression> {
    using _base = unary_expression<E, log_expression>;

public:
    constexpr explicit log_expression(const expression<E>& e) noexcept
        : _base{e()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        return std::log(_base::expr()(args...));
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        return _base::expr().template derivative<X>() / _base::expr();
    }
};

template<class E>
constexpr log_expression<E> log(const expression<E>& e) noexcept {
    return log_expression<E>{e()};
}

}

#endif
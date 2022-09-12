#ifndef SYMBOLIC_UNARY_EXPRESSION_HPP
#define SYMBOLIC_UNARY_EXPRESSION_HPP

#include "expression.hpp"

namespace SYMBOLIC_NAMESPACE {

template<class E, template<class, auto...> class Op, auto... Args>
class unary_expression : public expression<Op<E, Args...>> {
    const E _e;

public:
    constexpr explicit unary_expression(const expression<E>& e) noexcept
        : _e{e()} {}

    constexpr const E& expr() const noexcept {
        return _e;
    }
};

template<class E, template<class, auto...> class Op, auto... Args>
constexpr auto simplify(const unary_expression<E, Op, Args...>& e) noexcept {
    const auto se = simplify(e.expr());
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(e.expr())>, std::remove_cvref_t<decltype(se)>>)
        return Op<E, Args...>{simplify(e.expr())};
    else
        return simplify(Op<E, Args...>{simplify(e.expr())});
}

}

#endif
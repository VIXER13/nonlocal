#ifndef SYMBOLIC_BINARY_EXPRESSION_HPP
#define SYMBOLIC_BINARY_EXPRESSION_HPP

#include "expression.hpp"

#include <type_traits>
#include <utility>

namespace metamath::symbolic {

template<class E1, class E2, template<class, class> class Op>
class binary_expression : public expression<Op<E1, E2>> {
    const E1 _e1;
    const E2 _e2;

public:
    constexpr explicit binary_expression(const expression<E1>& e1, const expression<E2>& e2) noexcept
        : _e1{e1()}, _e2{e2()} {}

    constexpr std::pair<const E1&, const E2&> expr() const noexcept {
        return {_e1, _e2};
    }
};

template<class E1, class E2, template<class, class> class Op>
constexpr auto simplify(const binary_expression<E1, E2, Op>& e) noexcept {
    const auto [e1, e2] = e.expr();
    const auto [se1, se2] = std::pair{simplify(e1), simplify(e2)};
    if constexpr (std::is_same_v<std::remove_cvref_t<decltype(e1)>, std::remove_cvref_t<decltype(se1)>> &&
                  std::is_same_v<std::remove_cvref_t<decltype(e2)>, std::remove_cvref_t<decltype(se2)>>)
        return Op{se1, se2};
    else
        return simplify(Op{se1, se2});
}

}

#endif
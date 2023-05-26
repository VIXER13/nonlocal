#ifndef SYMBOLIC_SIMPLIFY_HPP
#define SYMBOLIC_SIMPLIFY_HPP

#include "expression.hpp"

#include <tuple>

namespace metamath::symbolic {

class _simplify final {
    constexpr explicit _simplify() noexcept = default;

    template<class Tuple, size_t... I>
    static constexpr auto simplify(const Tuple& expressions, const std::index_sequence<I...>) {
        return std::make_tuple(std::get<I>(expressions)...);
    }

public:
    template<class... E>
    friend constexpr auto simplify(const std::tuple<E...>& e);
};

template<class... E>
constexpr auto simplify(const std::tuple<E...>& e) {
    return _simplify::simplify(e, std::make_index_sequence<sizeof...(E)>{});
}

}

#endif
#pragma once

#include <metamath/symbolic/base/expression.hpp>

#include <cinttypes>
#include <tuple>
#include <utility>

namespace metamath::symbolic {

class _derivative final {
    constexpr explicit _derivative() noexcept = default;

    template<class E>
    static constexpr E derivative_impl(const E& e) {
        return e;
    }

    template<auto X, auto... Vars, class E>
    static constexpr auto derivative_impl(const E& e) {
        return derivative_impl<Vars...>(e.template derivative<X>());
    }

public:
    template<auto X, auto... Vars, class E>
    friend constexpr auto derivative(const E& e);

    template<auto X, auto... Vars, class... E>
    friend constexpr auto derivative(const std::tuple<E...>& e);
};

template<auto X, auto... Vars, class E>
constexpr auto derivative(const E& e) {
    return _derivative::derivative_impl<X, Vars...>(e);
}

template<auto X, auto... Vars, class... E>
constexpr auto derivative(const std::tuple<E...>& e) {
    return std::make_tuple(_derivative::derivative_impl<X, Vars...>(std::get<E>(e))...);
}

}
#pragma once

#include <metamath/symbolic/base/expression.hpp>

#include <tuple>

namespace metamath::symbolic {

class _basis_production final {
    constexpr explicit _basis_production() noexcept = default;

    template<class E, class... F, size_t... I>
    static constexpr auto product(const expression<E>& first, const std::tuple<F...>& second, const std::index_sequence<I...>) {
        return std::make_tuple(first() * std::get<I>(second)...);
    }

    template<class... E, class... F, size_t... I>
    static constexpr auto basis_production_impl(const std::tuple<E...>& first, const std::tuple<F...>& second, const std::index_sequence<I...>) {
        return std::tuple_cat(product(std::get<I>(first), second, std::make_index_sequence<sizeof...(F)>{})...);
    }

    template<class... E, class... F>
    static constexpr auto basis_production_impl(const std::tuple<E...>& first, const std::tuple<F...>& second) {
        return basis_production_impl(first, second, std::make_index_sequence<sizeof...(E)>{});
    }

    template<class... E>
    static constexpr std::tuple<E...> basis_production_impl(const std::tuple<E...> first) {
        return first;
    }

    template<class Basis, class... Bases>
    static constexpr auto basis_production_impl(const Basis& first, const Bases&... bases) {
        return basis_production_impl(first, basis_production_impl(bases...));
    }

public:
    template<class Basis, class... Bases>
    friend constexpr auto basis_production(const Basis& first, const Bases&... bases);
};

template<class Basis, class... Bases>
constexpr auto basis_production(const Basis& first, const Bases&... bases) {
    return _basis_production::basis_production_impl(first, _basis_production::basis_production_impl(bases...));
}

}
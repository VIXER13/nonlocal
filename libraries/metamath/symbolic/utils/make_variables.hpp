#pragma once

#include <metamath/symbolic/base/variable.hpp>

#include <tuple>

namespace metamath::symbolic {

class _make_variables final {
    constexpr explicit _make_variables() noexcept = default;

    template<size_t... I>
    static constexpr std::tuple<variable<I>...> make_variables_impl(const std::index_sequence<I...>&) noexcept {
        return std::make_tuple(variable<I>{}...);
    }

public:
    template<size_t N>
    friend constexpr auto make_variables() noexcept;
};

template<size_t N>
constexpr auto make_variables() noexcept {
    return _make_variables::make_variables_impl(std::make_index_sequence<N>{});
}

}
#ifndef SYMDIFF_MAKE_VARIABLES_HPP
#define SYMDIFF_MAKE_VARIABLES_HPP

#include "variable.hpp"
#include <tuple>

namespace metamath::symdiff {

class _make_variables final {
    constexpr explicit _make_variables() noexcept = default;

    template<uintmax_t... I>
    static constexpr std::tuple<variable<I>...> make_variables_impl(const std::index_sequence<I...>&) noexcept {
        return std::make_tuple(variable<I>{}...);
    }

public:
    template<uintmax_t N>
    friend constexpr auto make_variables() noexcept;
};

template<uintmax_t N>
constexpr auto make_variables() noexcept {
    return _make_variables::make_variables_impl(std::make_index_sequence<N>{});
}

}

#endif
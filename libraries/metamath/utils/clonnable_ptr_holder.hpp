#pragma once

#include <concepts>
#include <memory>
#include <tuple>

namespace metamath::utils {

template<class T>
concept clonnable = requires(const T& v) { { v.clone() } -> std::convertible_to<std::unique_ptr<T>>; };

template<class... Ts>
requires(sizeof...(Ts) > 0) && (clonnable<Ts> && ...)
class clonnable_ptrs {
    using tuple_type = std::tuple<std::unique_ptr<Ts>...>;

    static tuple_type _clone(const tuple_type& src) {
        constexpr auto clone = [](const auto&... ptrs) { return tuple_type{(ptrs ? ptrs->clone() : nullptr)...}; };
        return std::apply(clone, src);
    }

public:
    tuple_type _ptrs{};

public:
    constexpr clonnable_ptrs() noexcept = default;

    constexpr explicit clonnable_ptrs(std::unique_ptr<Ts>... ptrs) noexcept
        : _ptrs(std::move(ptrs)...) {}

    clonnable_ptrs(const clonnable_ptrs& other)
        : _ptrs(_clone(other._ptrs)) {}
};

} // namespace metamath::utils

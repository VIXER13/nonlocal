#pragma once

#include <concepts>
#include <memory>
#include <tuple>

namespace metamath::types {

template<class T>
concept copyable = requires(const T& v) { { v.copy() } -> std::convertible_to<std::unique_ptr<T>>; };

template<class... Ts>
requires(sizeof...(Ts) > 0) && (copyable<Ts> && ...)
class copyable_ptrs {
    using tuple_type = std::tuple<std::unique_ptr<Ts>...>;

    static tuple_type _copy(const tuple_type& src) {
        constexpr auto copy = [](const auto&... ptrs) { return tuple_type{(ptrs ? ptrs->copy() : nullptr)...}; };
        return std::apply(copy, src);
    }

public:
    tuple_type _ptrs{};

public:
    constexpr copyable_ptrs() noexcept = default;

    constexpr explicit copyable_ptrs(std::unique_ptr<Ts>... ptrs) noexcept
        : _ptrs(std::move(ptrs)...) {}

    copyable_ptrs(const copyable_ptrs& other)
        : _ptrs(_copy(other._ptrs)) {}
};

} // namespace metamath::types

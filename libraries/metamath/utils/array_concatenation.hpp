#ifndef METAMATH_ARRAY_CONCATENATION_HPP
#define METAMATH_ARRAY_CONCATENATION_HPP

#include <array>
#include <cstddef>
#include <utility>

namespace metamath::utils {

class _concat final {
    constexpr explicit _concat() noexcept = default;

    template<class T, size_t N, size_t M, size_t... I, size_t... J>
    static constexpr std::array<T, N + M> concat_impl(const std::array<T, N>& first,   const std::array<T, M>& second,
                                                 const std::index_sequence<I...>, const std::index_sequence<J...>) {
        return {first[I]..., second[J]...};
    }

    template<class T, size_t N, size_t M>
    static constexpr std::array<T, N + M> concat_impl(const std::array<T, N>& first, const std::array<T, M>& second) {
        return concat_impl(first, second, std::make_index_sequence<N>{}, std::make_index_sequence<M>{});
    }

    template<class T, size_t N, size_t... I>
    static constexpr std::array<T, N + (I + ...)> concat_impl(const std::array<T, N>& first, const std::array<T, I>&... arrays) {
        return concat_impl(first, concat_impl(arrays...));
    }

public:
    template<class T, size_t N, size_t... I>
    friend constexpr std::array<T, N + (I + ...)> concat(const std::array<T, N>& first, const std::array<T, I>&... arrays);
};

template<class T, size_t N, size_t... I>
constexpr std::array<T, N + (I + ...)> concat(const std::array<T, N>& first, const std::array<T, I>&... arrays) {
    return _concat::concat_impl(first, arrays...);
}

template<class T, size_t N>
constexpr std::array<T, N> concat(const std::array<T, N> first) {
    return first;
}

}

#endif
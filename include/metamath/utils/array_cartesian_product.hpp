#ifndef METAMATH_ARRAY_CARTESIAN_PRODUCT_HPP
#define METAMATH_ARRAY_CARTESIAN_PRODUCT_HPP

#include "array_concatenation.hpp"

namespace metamath::utils {

class _array_cartesian_product final {
    constexpr explicit _array_cartesian_product() noexcept = default;

    template<class T, size_t M, size_t... I>
    static constexpr std::array<std::array<T, 2>, M>
    generate_arrays(const T val, const std::array<T, M>& arr, const std::index_sequence<I...>) {
        return {std::array{val, arr[I]}...};
    }

    template<class T, size_t K, size_t M, size_t... I>
    static constexpr std::array<std::array<T, K + 1>, M>
    generate_arrays(const T val, const std::array<std::array<T, K>, M>& arr, const std::index_sequence<I...>) {
        return std::array{concat(std::array{val}, arr[I])...};
    }

    template<class T, size_t N, size_t M, size_t... I>
    static constexpr std::array<std::array<T, 2>, N * M>
    cartesian_product(const std::array<T, N>& first, const std::array<T, M>& second, const std::index_sequence<I...>) {
        return {concat(generate_arrays(first[I], second, std::make_index_sequence<M>{})...)};
    }

    template<class T, size_t N, size_t K, size_t M, size_t... I>
    static constexpr std::array<std::array<T, K + 1>, N * M>
    cartesian_product(const std::array<T, N>& first, const std::array<std::array<T, K>, M>& second, const std::index_sequence<I...>) {
        return {concat(generate_arrays(first[I], second, std::make_index_sequence<M>{})...)};
    }

    template<class T, size_t N, size_t M>
    static constexpr std::array<std::array<T, 2>, N * M>
    cartesian_product(const std::array<T, N>& first, const std::array<T, M>& second) {
        return cartesian_product(first, second, std::make_index_sequence<N>{});
    }

    template<class T, size_t N, size_t K, size_t M>
    static constexpr std::array<std::array<T, K + 1>, N * M>
    cartesian_product(const std::array<T, N>& first, const std::array<std::array<T, K>, M>& second) {
        return cartesian_product(first, second, std::make_index_sequence<N>{});
    }

    template<class T, size_t N, size_t... I>
    static constexpr std::array<std::array<T, sizeof...(I) + 1>, N * (I * ...)>
    cartesian_product_impl(const std::array<T, N>& first, const std::array<T, I>&... arrays) {
        return cartesian_product(first, cartesian_product(arrays...));
    }

public:
    template<class T, size_t N, size_t... I>
    friend constexpr std::array<std::array<T, sizeof...(I) + 1>, N * (I * ...)>
    array_cartesian_product(const std::array<T, N>& first, const std::array<T, I>&... arrays);
};

template<class T, size_t N, size_t... I>
constexpr std::array<std::array<T, sizeof...(I) + 1>, N * (I * ...)>
array_cartesian_product(const std::array<T, N>& first, const std::array<T, I>&... arrays) {
    return _array_cartesian_product::cartesian_product(first, arrays...);
}

template<class T, size_t N>
constexpr std::array<T, N> cartesian_product(const std::array<T, N> first) {
    return first;
}

}

#endif
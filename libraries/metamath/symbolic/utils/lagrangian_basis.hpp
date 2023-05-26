#ifndef SYMBOLIC_LAGRANGIAN_BASIS_HPP
#define SYMBOLIC_LAGRANGIAN_BASIS_HPP

#include "variable.hpp"

namespace metamath::symbolic {

class _lagrangian_function final {
    explicit constexpr _lagrangian_function() noexcept = default;

    template<size_t X, size_t K, size_t I, class T, size_t N>
    static constexpr auto generate_term(const std::array<T, N>& nodes) noexcept {
        if constexpr(K == I)
            return metamath::symbolic::integral_constant<1>{};
        else {
            constexpr variable<X> x;
            return (x - nodes[I]) / (nodes[K] - nodes[I]);
        }
    }

    template<size_t X, size_t K, class T, size_t N, size_t... I>
    static constexpr auto generate_lagrangian_function(const std::array<T, N>& nodes, const std::index_sequence<I...>) noexcept {
        return (generate_term<X, K, I>(nodes) * ...);
    }

public:
    template<size_t X, size_t K, class T, size_t N>
    friend constexpr auto generate_lagrangian_function(const std::array<T, N>& nodes) noexcept;
};

template<size_t X, size_t K, class T, size_t N>
constexpr auto generate_lagrangian_function(const std::array<T, N>& nodes) noexcept {
    static_assert(N > 0, "The number of nodes must be greater than 0");
    return _lagrangian_function::generate_lagrangian_function<X, K>(nodes, std::make_index_sequence<N>{});
}

class _lagrangian_basis final {
    explicit constexpr _lagrangian_basis() noexcept = default;

    template<size_t X, class T, size_t N, size_t... K>
    static constexpr auto generate_lagrangian_basis(const std::array<T, N>& nodes, const std::index_sequence<K...>) noexcept {
        return std::make_tuple(generate_lagrangian_function<X, K>(nodes)...);
    }

public:
    template<size_t X, class T, size_t N>
    friend constexpr auto generate_lagrangian_basis(const std::array<T, N>& nodes) noexcept;
};

template<size_t X, class T, size_t N>
constexpr auto generate_lagrangian_basis(const std::array<T, N>& nodes) noexcept {
    return _lagrangian_basis::generate_lagrangian_basis<X>(nodes, std::make_index_sequence<N>{});
}

}

#endif
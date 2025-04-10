#pragma once

#include "plus.hpp"

namespace metamath::symbolic {

template<class E1, class E2>
class multiplies : public binary_expression<E1, E2, multiplies> {
    using _base = binary_expression<E1, E2, multiplies>;

public:
    constexpr explicit multiplies(const expression<E1>& e1, const expression<E2>& e2) noexcept
        : _base{e1(), e2()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        const auto [e1, e2] = _base::expr();
        return e1(args...) * e2(args...);
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        const auto [e1, e2] = _base::expr();
        return e1.template derivative<X>() * e2 + e1 * e2.template derivative<X>();
    }
};

template<class E1, class E2>
constexpr multiplies<E1, E2> operator*(const expression<E1>& e1, const expression<E2>& e2) noexcept {
    return multiplies<E1, E2>{e1(), e2()};
}

template<class T, class E, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
constexpr multiplies<constant<T>, E> operator*(const T& c, const expression<E>& e) noexcept {
    return multiplies<constant<T>, E>{constant<T>{c}, e()};
}

template<class E, class T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
constexpr multiplies<E, constant<T>> operator*(const expression<E>& e, const T& c) noexcept {
    return multiplies<E, constant<T>>{e(), constant<T>{c}};
}

template<auto N1, auto N2>
constexpr integral_constant<N1 * N2> simplify(const multiplies<integral_constant<N1>, integral_constant<N2>>) noexcept {
    return {};
}

template<auto N, class T>
constexpr constant<decltype(N * T{})> simplify(const multiplies<integral_constant<N>, constant<T>>& e) noexcept {
    const auto [_, c] = e.expr();
    return {N * c()};
}

template<class T, auto N>
constexpr constant<decltype(T{} * N)> simplify(const multiplies<constant<T>, integral_constant<N>>& e) noexcept {
    const auto [c, _] = e.expr();
    return {c() * N};
}

template<class T1, class T2>
constexpr constant<decltype(T1{} + T2{})> simplify(const multiplies<constant<T1>, constant<T2>>& e) noexcept {
    const auto [c1, c2] = e.expr();
    return {c1() + c2()};
}

template<auto N, class E, std::enable_if_t<N == 0, bool> = true>
constexpr integral_constant<N> simplify(const multiplies<integral_constant<N>, E>&) noexcept {
    return {};
}

template<class E, auto N, std::enable_if_t<N == 0, bool> = true>
constexpr integral_constant<N> simplify(const multiplies<E, integral_constant<N>>&) noexcept {
    return {};
}

template<auto N, class E, std::enable_if_t<N == 1, bool> = true>
constexpr auto simplify(const multiplies<integral_constant<N>, E>& m) noexcept {
    const auto [_, e] = m.expr();
    return simplify(e);
}

template<class E, auto N, std::enable_if_t<N == 1, bool> = true>
constexpr auto simplify(const multiplies<E, integral_constant<N>>& m) noexcept {
    const auto [e, _] = m.expr();
    return simplify(e);
}

}
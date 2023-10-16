#ifndef SYMBOLIC_DIVIDES_HPP
#define SYMBOLIC_DIVIDES_HPP

#include "minus.hpp"
#include "multiplies.hpp"
#include "power_expression.hpp"

namespace metamath::symbolic {

template<class E1, class E2>
class divides : public binary_expression<E1, E2, divides> {
    using _base = binary_expression<E1, E2, divides>;

public:
    constexpr explicit divides(const expression<E1>& e1, const expression<E2>& e2) noexcept
        : _base{e1(), e2()} {}

    template<class... Args>
    constexpr auto operator()(const Args&... args) const {
        const auto [e1, e2] = _base::expr();
        return e1(args...) / e2(args...);
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        const auto [e1, e2] = _base::expr();
        return (e1.template derivative<X>() * e2 - e1 * e2.template derivative<X>()) / power<2>(e2);
    }
};

template<class E1, class E2>
constexpr divides<E1, E2> operator/(const expression<E1>& e1, const expression<E2>& e2) noexcept {
    return divides<E1, E2>{e1(), e2()};
}

template<class T, class E, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
constexpr divides<constant<T>, E> operator/(const T& c, const expression<E>& e) noexcept {
    return divides<constant<T>, E>{constant<T>{c}, e()};
}

template<class E, class T, std::enable_if_t<std::is_arithmetic_v<T>, bool> = true>
constexpr divides<E, constant<T>> operator/(const expression<E>& e, const T& c) noexcept {
    return divides<E, constant<T>>{e(), constant<T>{c}};
}

template<auto N1, auto N2>
constexpr integral_constant<N1 / N2> simplify(const divides<integral_constant<N1>, integral_constant<N2>>) noexcept {
    return {};
}

template<auto N, class T>
constexpr constant<decltype(N / T{})> simplify(const divides<integral_constant<N>, constant<T>>& e) noexcept {
    const auto [_, c] = e.expr();
    return {N / c()};
}

template<class T, auto N>
constexpr constant<decltype(T{} / N)> simplify(const divides<constant<T>, integral_constant<N>>& e) noexcept {
    const auto [c, _] = e.expr();
    return {c() / N};
}

template<auto N, class E, std::enable_if_t<N == 0, bool> = true>
constexpr integral_constant<N> simplify(const divides<integral_constant<N>, E>&) noexcept {
    return {};
}

template<class E, auto N, std::enable_if_t<N == 1 || N == 0, bool> = true>
constexpr auto simplify(const divides<E, integral_constant<N>>& d) noexcept {
    static_assert(N != 0, "Division by zero.");
    const auto [e, _] = d.expr();
    return simplify(e);
}

template<class E1U, class E1D, class E2U, class E2D>
constexpr auto simplify(const multiplies<divides<E1U, E1D>, divides<E2U, E2D>>& e) noexcept {
    const auto [l, r] = e.expr();
    const auto [lu, ld] = l.expr();
    const auto [ru, rd] = r.expr();
    return simplify((lu * ru) / (ld * rd));
}

}

#endif
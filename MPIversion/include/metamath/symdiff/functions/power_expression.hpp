#ifndef SYMDIFF_POWER_EXPRESSION_HPP
#define SYMDIFF_POWER_EXPRESSION_HPP

#include "power.hpp"
#include "multiplies.hpp"

namespace metamath::symdiff {

template<class E, intmax_t N>
struct power_expression;

template<intmax_t N, class E>
constexpr std::enable_if_t<!N, integral_constant<intmax_t, 1>> power(const expression<E>& e) {
    return integral_constant<intmax_t, 1>{};
}

template<intmax_t N, class E>
constexpr std::enable_if_t<N == 1, const E&> power(const expression<E>& e) {
    return e();
}

template<intmax_t N, class E>
constexpr std::enable_if_t<N && N != 1, power_expression<E, N>> power(const expression<E>& e);

template<intmax_t N, class E, intmax_t M>
constexpr std::enable_if_t<N && N != 1, power_expression<E, N * M>> power(const power_expression<E, M>& e);

template<typename E>
struct is_pow_expression : std::false_type {};

template<typename E, intmax_t N>
struct is_pow_expression<power_expression<E, N>> : std::true_type {};

template<class E>
struct pow_expression_index : std::integral_constant<intmax_t, 0> {};

template<class E, intmax_t N>
struct pow_expression_index<power_expression<E, N>> : std::integral_constant<intmax_t, N> {};

template<class E, intmax_t N>
using power_expression_type = std::conditional_t<
    !N,
    integral_constant<intmax_t, 1>,
    std::conditional_t<
        N == 1,
        E,
        std::conditional_t<
            is_pow_expression<E>{},
            power_expression<E, N * pow_expression_index<E>{}>,
            power_expression<E, N>
        >
    >
>;

template<class E, intmax_t N>
class power_expression : public expression<power_expression<E, N>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = multiplies_type<
        multiplies_type<integral_constant<intmax_t, N>, power_expression_type<E, N-1>>,
        typename E::template derivative_type<X>
    >;

    constexpr explicit power_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(function::power<N>(e(x))) {
        return function::power<N>(e(x));
    }

    constexpr const E& expression() const { return e; }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return integral_constant<intmax_t, N>{} * power<N-1>(e) * e.template derivative<X>();
    }
};

template<intmax_t N, class E>
constexpr std::enable_if_t<N && N != 1, power_expression<E, N>> power(const expression<E>& e) {
    return power_expression<E, N>{e};
}

template<intmax_t N, class E, intmax_t M>
constexpr std::enable_if_t<N && N != 1, power_expression<E, N * M>> power(const power_expression<E, M>& e) {
    return power_expression<E, N * M>{e.expression()};
}

}

#endif
#ifndef METAMATH_CONSTANTS_MULTINOMIAL_HPP
#define METAMATH_CONSTANTS_MULTINOMIAL_HPP

#include "factorial.hpp"

namespace metamath::constants {

// Мультиномиальные коэффициенты
class _multinomial {
    template<class Type, uintmax_t Value, uintmax_t Current_K, uintmax_t... K>
    struct calc_multinomial : calc_multinomial<Type, Value / factorial<Current_K, Type>{}, K...> {};

    template<class Type, uintmax_t Value, uintmax_t K>
    struct calc_multinomial<Type, Value, K> : std::integral_constant<Type, Value / factorial<K, Type>{}> {};

public:
    template<uintmax_t... K>
    friend struct multinomial;
};

template<uintmax_t... K>
struct multinomial : _multinomial::calc_multinomial<intmax_t, factorial<(K + ...), intmax_t>{}, K...> {};

template<>
struct multinomial<> : std::integral_constant<intmax_t, 1> {};

// Биномиальные коэффициенты хоть и являются частным случаем мультиномиальных,
// но такое их вычисление позволяет вычислять коэффициенты для больших степеней биномов.
template<uintmax_t N, uintmax_t K, class Type = intmax_t>
struct binomial : std::integral_constant<
    std::enable_if_t<
            binomial<N-1, K, Type>{} < std::numeric_limits<Type>::max() - binomial<N-1, K-1, Type>{},
            Type
        >,
        binomial<N-1, K-1, Type>{} + binomial<N-1, K, Type>{}
    > {};

template<uintmax_t K, class Type>
struct binomial<0, K, Type> : std::integral_constant<Type, 0> {};

template<class Type>
struct binomial<0, 0, Type> : std::integral_constant<Type, 1> {};

}

#endif
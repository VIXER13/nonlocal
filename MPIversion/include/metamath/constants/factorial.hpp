#ifndef METAMATH_CONSTANTS_FACTORIAL_HPP
#define METAMATH_CONSTANTS_FACTORIAL_HPP

#include "prime_number.hpp"

namespace metamath::constants {

// M-кратный неполный факториал. Вычисляется как произведение чисел в интервале [N, L) с шагом -M.
template<uintmax_t N, uintmax_t L, uintmax_t M, class Type = intmax_t, class Enable = void>
struct incomplete_multifactorial;

template<uintmax_t N, uintmax_t L, uintmax_t M, class Type>
struct incomplete_multifactorial<N, L, M, Type, std::enable_if_t<N < M || N <= L>> :
    std::integral_constant<Type, 1> {};

template<uintmax_t N, uintmax_t L, uintmax_t M, class Type>
struct incomplete_multifactorial<N, L, M, Type, std::enable_if_t<(N >= M) && (N > L)>> :
    std::integral_constant<
        std::enable_if_t<
            (std::numeric_limits<Type>::max() / incomplete_multifactorial<N-M, L, M, Type>{} / N > 0),
            Type
        >,
        N * incomplete_multifactorial<N-M, L, M, Type>{}
    > {};

template<uintmax_t N, uintmax_t M, class Type = intmax_t>
using multifactorial = incomplete_multifactorial<N, 0, M, Type>;

template<uintmax_t N, uintmax_t L, class Type = intmax_t>
using incomplete_factorial = incomplete_multifactorial<N, L, 1, Type>;

template<uintmax_t N, class Type = intmax_t>
using factorial = incomplete_multifactorial<N, 0, 1, Type>;

// Субфакториал или количество беспорядков порядка N. В математике обозначают !N.
class _subfactorial {
    _subfactorial() = delete;

    template<uintmax_t K, uintmax_t N, class Type>
    struct calc_subfactorial :
        std::integral_constant<
            Type,
            (K % 2 ? -1 : 1) * incomplete_factorial<N, K, Type>{} + calc_subfactorial<K-1, N, Type>{}
        > {};

    template<uintmax_t N, class Type>
    struct calc_subfactorial<0, N, Type> : std::integral_constant<Type, factorial<N, Type>{}> {};

public:
    template<uintmax_t N, class Type>
    friend struct subfactorial;
};

template<uintmax_t N, class Type = intmax_t>
struct subfactorial : _subfactorial::calc_subfactorial<N, N, Type> {};

template<class Type>
struct subfactorial<0, Type> : std::integral_constant<Type, 1> {};

// Праймориал - произведение первых N простых чисел.
template<uintmax_t N, class Type = uintmax_t>
struct primorial :
    std::integral_constant<
        std::enable_if_t<
            (std::numeric_limits<Type>::max() / prime_at<N-1, Type>{} / primorial<N-1, Type>{} > 0),
            Type
        >,
        prime_at<N-1, Type>{} * primorial<N-1, Type>{}
    > {};

template<class Type>
struct primorial<1, Type> : std::integral_constant<Type, prime_at<0, Type>{}> {};

template<class Type>
struct primorial<0, Type> : std::integral_constant<Type, 1> {};

}

#endif
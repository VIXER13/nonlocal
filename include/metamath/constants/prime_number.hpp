#ifndef METAMATHTEST_PRIME_NUMBER_HPP
#define METAMATHTEST_PRIME_NUMBER_HPP

#include <limits>
#include <type_traits>

namespace metamath::constants {

class _prime_numbers {
    _prime_numbers() = delete;

    template<uintmax_t N, uintmax_t M, class Type>
    struct is_prime_check : std::conditional_t<
        N % M == 0,
        std::false_type,
        is_prime_check<N, M - 1, Type>
    > {};

    template<uintmax_t N, class Type>
    struct is_prime_check<N, 1, Type> : std::true_type {};

    template<uintmax_t N, class Type>
    struct find : std::conditional_t<
        is_prime_check<N, N / 2, Type>{},
        std::integral_constant<Type, N>,
        std::enable_if_t<(std::numeric_limits<Type>::max() - N > 1), find<N+2, Type>>
    > {};

public:
    // Является ли число N простым.
    template<uintmax_t N, class Type>
    friend struct is_prime;

    // Поиск простого числа под номером N (нумерация с 0).
    template<uintmax_t N, class Type>
    friend struct prime_at;
};

template<uintmax_t N, class Type = uintmax_t>
struct is_prime : _prime_numbers::is_prime_check<N, N / 2, Type> {};

template<uintmax_t N, class Type = uintmax_t>
struct prime_at : _prime_numbers::find<
    prime_at<N - 1, Type>{} + 2,
    std::enable_if_t<(std::numeric_limits<Type>::max() - prime_at<N - 1, Type>{} > 1), Type>
> {};

template<class Type>
struct prime_at<0, Type> : std::integral_constant<Type, 2> {};

template<class Type>
struct prime_at<1, Type> : std::integral_constant<Type, 3> {};

}

#endif
#ifndef METAMATH_FUNCTIONS_HERMITE_HPP
#define METAMATH_FUNCTIONS_HERMITE_HPP

#include <utility>
#include <stdexcept>
#include "power.hpp"
#include "factorial.hpp"

namespace metamath::function {

class _hermite {
    _hermite() = delete;

    template<uintmax_t N, class T>
    static constexpr std::enable_if_t<N == 0, std::pair<T, T>> hermite_term(const T x) {
        return {0, 1};
    }

    template<uintmax_t N, class T>
    static constexpr std::enable_if_t<N == 1, std::pair<T, T>> hermite_term(const T x) {
        return {1, x + x};
    }

    template<uintmax_t N, class T>
    static constexpr std::enable_if_t<(N > 1), std::pair<T, T>> hermite_term(const T x) {
        const auto [H2, H1] = hermite_term<N-1>(x);
        return {H1, (x + x) * H1 - (2 * N - 2) * H2};
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<K == 0, T> hermite_fast_term(const T x) {
        return power<N>(x);
    }

    template<uintmax_t N, uintmax_t K, class T>
    static constexpr std::enable_if_t<(K > 0), T> hermite_fast_term(const T x) {
        static_assert(!std::is_unsigned_v<T>, "The T must be signed integer or floating point.");
        constexpr int64_t coeff = (K % 2 ? -1 : 1) * constants::incomplete_factorial<N, K>{} / constants::factorial<N-2*K>{};
        return coeff * power<N-2*K>(x) + hermite_fast_term<N, K-1>(x);
    }

public:
    template<uintmax_t N, class T>
    friend constexpr T hermite(const T x);

    template<uintmax_t N, class T, bool With_Check>
    friend constexpr T hermite_fast(const T x);
};

    template<uintmax_t N, class T>
    constexpr T hermite(const T x) {
        return _hermite::hermite_term<N>(x).second;
    }

    template<uintmax_t N, class T, bool With_Check>
    constexpr T hermite_fast(const T x) {
        if constexpr (With_Check) {
            if(x < -1 || x > 1)
                throw std::domain_error{"Argument out of range in hermite_fast."};
        }
        return _hermite::hermite_fast_term<N, N/2>(x + x);
    }

}

#endif
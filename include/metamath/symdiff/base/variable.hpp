#ifndef SYMDIFF_VARIABLE_HPP
#define SYMDIFF_VARIABLE_HPP

#include "integral_constant.hpp"

namespace metamath::symdiff {

template<uintmax_t N>
struct variable : expression<variable<N>> {
    template<uintmax_t X>
    using derivative_type = integral_constant<intmax_t, N == X>;

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(x[N]) {
        return x[N];
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return derivative_type<X>{};
    }

    constexpr operator uintmax_t() const {
        return N;
    }
};

}

#endif
#ifndef SYMDIFF_EXP_EXPRESSION_HPP
#define SYMDIFF_EXP_EXPRESSION_HPP

#include <cmath>
#include "multiplies.hpp"

namespace metamath::symdiff {

template<class E>
class exp_expression : public expression<exp_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = multiplies_type<exp_expression<E>, typename E::template derivative_type<X>>;

    constexpr explicit exp_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::exp(e(x))) {
        return std::exp(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return exp(e) * e.template derivative<X>();
    }
};

template<class E>
constexpr exp_expression<E> exp(const expression<E>& e) {
    return exp_expression<E>{e};
}

}

#endif
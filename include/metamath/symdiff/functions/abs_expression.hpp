#ifndef METAMATHTEST_ABS_EXPRESSION_HPP
#define METAMATHTEST_ABS_EXPRESSION_HPP

#include <cmath>
#include "multiplies.hpp"
#include "sign_expression.hpp"

namespace metamath::symdiff {

template<class E>
class abs_expression : public expression<abs_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = multiplies_type<
        typename E::template derivative_type<X>,
        sign_expression<E>
    >;

    constexpr abs_expression(const expression<E> &e) : e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::abs(e(x))) {
        return std::abs(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return e.template derivative<X>() * sign(e);
    }
};

template<class E>
constexpr abs_expression<E> abs(const expression<E>& e) {
    return abs_expression<E>{e};
}

}

#endif
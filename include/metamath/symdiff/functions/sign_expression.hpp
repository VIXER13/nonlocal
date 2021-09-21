#ifndef SYMDIFF_SIGN_EXPRESSION_HPP
#define SYMDIFF_SIGN_EXPRESSION_HPP

#include "integral_constant.hpp"

namespace metamath::symdiff {

template<class E>
class sign_expression : public expression<sign_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = integral_constant<intmax_t, 0>;

    constexpr explicit sign_expression(const expression<E>& e) : e{e()} {}

    template<class U>
    constexpr intmax_t operator()(const U& x) const {
        const auto value = e(x);
        return value < 0 ? -1 :
               value > 0 ?  1 : 0;
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return derivative_type<X>{};
    }
};

template<class E>
constexpr sign_expression<E> sign(const expression<E>& e) {
    return sign_expression<E>{e};
}

}

#endif
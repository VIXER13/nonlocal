#ifndef SYMDIFF_TAN_EXPRESSION_HPP
#define SYMDIFF_TAN_EXPRESSION_HPP

#include <cmath>
#include "divides.hpp"
#include "cos_expression.hpp"
#include "power_expression.hpp"

namespace metamath::symdiff {

template<class E>
class tan_expression : public expression<tan_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = divides_type<
        typename E::template derivative_type<X>,
        power_expression_type<cos_expression<E>, 2>
    >;

    constexpr explicit tan_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::tan(e(x))) {
        return std::tan(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return e.template derivative<X>() / power<2>(cos(e));
    }
};

template<class E>
constexpr tan_expression<E> tan(const expression<E>& e) {
    return tan_expression<E>{e};
}

}

#endif
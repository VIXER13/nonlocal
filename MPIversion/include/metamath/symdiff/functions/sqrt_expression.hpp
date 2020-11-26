#ifndef SYMDIFF_SQRT_EXPRESSION_HPP
#define SYMDIFF_SQRT_EXPRESSION_HPP

#include <cmath>
#include "divides.hpp"

namespace metamath::symdiff {

template<class E>
class sqrt_expression : public expression<sqrt_expression<E>> {
    const E e;

    template<uintmax_t X>
    struct derivative_t {
        using derivative_type = divides_type<
            typename E::template derivative_type<X>,
            multiplies_type<integral_constant<int, 2>, sqrt_expression<E>>
        >;
    };

public:
    template<uintmax_t X>
    using derivative_type = typename derivative_t<X>::derivative_type;

    constexpr explicit sqrt_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::sqrt(e(x))) {
        return std::sqrt(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return e.template derivative<X>() / (integral_constant<int, 2>{} * sqrt(e));
    }
};

template<class E>
constexpr sqrt_expression<E> sqrt(const expression<E>& e) {
    return sqrt_expression<E>{e};
}

}

#endif
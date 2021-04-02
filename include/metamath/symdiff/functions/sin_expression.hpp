#ifndef SYMDIFF_SIN_EXPRESSION_HPP
#define SYMDIFF_SIN_EXPRESSION_HPP

#include <cmath>
#include "multiplies.hpp"

namespace metamath::symdiff {

template<class E>
struct cos_expression;

template<class E>
class sin_expression : public expression<sin_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = multiplies_type<cos_expression<E>, typename E::template derivative_type<X>>;

    constexpr explicit sin_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::sin(e(x))) {
        return std::sin(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return cos(e) * e.template derivative<X>();
    }
};

template<class E>
constexpr sin_expression<E> sin(const expression<E>& e){
    return sin_expression<E>{e};
}

}

#endif
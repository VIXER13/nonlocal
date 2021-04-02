#ifndef SYMDIFF_COS_EXPRESSION_HPP
#define SYMDIFF_COS_EXPRESSION_HPP

#include <cmath>
#include "multiplies.hpp"
#include "negate.hpp"

namespace metamath::symdiff {

template<class E>
struct sin_expression;

template<class E>
class cos_expression : public expression<cos_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = negate_type<multiplies_type<sin_expression<E>, typename E::template derivative_type<X>>>;

    constexpr explicit cos_expression(const expression<E>& e) :
        e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::cos(e(x))) {
        return std::cos(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return -(sin(e) * e.template derivative<X>());
    }
};

template<class E>
constexpr cos_expression<E> cos(const expression<E>& e){
    return cos_expression<E>{e};
}

}

#endif
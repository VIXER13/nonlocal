#ifndef SYMDIFF_LOG_EXPRESSION_HPP
#define SYMDIFF_LOG_EXPRESSION_HPP

#include <cmath>
#include "divides.hpp"

namespace metamath::symdiff {

template<class E>
class log_expression : public expression<log_expression<E>> {
    const E e;

public:
    template<uintmax_t X>
    using derivative_type = divides_type<typename E::template derivative_type<X>, E>;

    constexpr explicit log_expression(const expression<E>& e) : e{e()} {}

    template<class U>
    constexpr auto operator()(const U& x) const -> decltype(std::log(e(x))) {
        return std::log(e(x));
    }

    template<uintmax_t X>
    constexpr derivative_type<X> derivative() const {
        return e.template derivative<X>() / e;
    }
};

template<class E>
constexpr log_expression<E> log(const expression<E>& e) {
    return log_expression<E>{e};
}

}

#endif
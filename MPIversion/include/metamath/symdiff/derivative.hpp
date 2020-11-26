#ifndef SYMDIFF_DERIVATIVE_HPP
#define SYMDIFF_DERIVATIVE_HPP

#include <tuple>
#include <utility>

namespace metamath::symdiff {

template<uintmax_t X, class E>
constexpr typename E::template derivative_type<X> derivative(const E& e) {
    return e.template derivative<X>();
}

class _derivative_tuple {
    constexpr explicit _derivative_tuple() noexcept = default;

    template<uintmax_t X, class Tuple, size_t... I>
    static constexpr auto derivative_tuple(const Tuple& expressions, const std::index_sequence<I...>&) {
        return std::make_tuple(derivative<X>(std::get<I>(expressions))...);
    }

public:
    template<uintmax_t X, class... E>
    friend constexpr std::tuple<typename E::template derivative_type<X>...> derivative(const std::tuple<E...>& e);
};

template<uintmax_t X, class... E>
constexpr std::tuple<typename E::template derivative_type<X>...> derivative(const std::tuple<E...>& e) {
    return _derivative_tuple::derivative_tuple<X>(e, std::make_index_sequence<sizeof...(E)>{});
}

}

#endif
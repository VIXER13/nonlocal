#ifndef SYMDIFF_DERIVATIVE_HPP
#define SYMDIFF_DERIVATIVE_HPP

#include <tuple>
#include <utility>

namespace metamath::symdiff {

class _derivative {
    constexpr explicit _derivative() noexcept = default;

    template<class E>
    static constexpr E derivative_impl(const E& e) {
        return e;
    }

    template<uintmax_t X, uintmax_t... Vars, class E>
    static constexpr auto derivative_impl(const E& e) {
        return derivative_impl<Vars...>(e.template derivative<X>());
    }

    template<uintmax_t... Vars, class Tuple, size_t... I>
    static constexpr auto derivative_tuple_impl(const Tuple& expressions, const std::index_sequence<I...>&) {
        return std::make_tuple(derivative_impl<Vars...>(std::get<I>(expressions))...);
    }

public:
    template<uintmax_t X, uintmax_t... Vars, class E>
    friend constexpr auto derivative(const E& e);

    template<uintmax_t X, uintmax_t... Vars, class... E>
    friend constexpr auto derivative(const std::tuple<E...>& e);
};

template<uintmax_t X, uintmax_t... Vars, class E>
constexpr auto derivative(const E& e) {
    return _derivative::derivative_impl<X, Vars...>(e);
}

template<uintmax_t X, uintmax_t... Vars, class... E>
constexpr auto derivative(const std::tuple<E...>& e) {
    return _derivative::derivative_tuple_impl<X, Vars...>(e, std::make_index_sequence<sizeof...(E)>{});
}

}

#endif
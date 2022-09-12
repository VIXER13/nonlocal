#ifndef SYMDIFF_TO_FUNCTION_HPP
#define SYMDIFF_TO_FUNCTION_HPP

#include <array>
#include <tuple>
#include <functional>

namespace SYMBOLIC_NAMESPACE {

template<class T, size_t N, class E>
std::function<T(const std::array<T, N>&)> to_function(const E& e) {
    return [e](const std::array<T, N>& x) { return e(x); };
}

class _to_array_of_functions final {
    constexpr explicit _to_array_of_functions() noexcept = default;

    template<class T, size_t N, class Tuple, size_t... I>
    static std::array<std::function<T(const std::array<T, N>&)>, sizeof...(I)>
    to_array_of_functions(const Tuple& expressions, const std::index_sequence<I...>) {
        return {to_function<T, N>(std::get<I>(expressions))...};
    }

public:
    template<class T, size_t N, class... E>
    friend std::array<std::function<T(const std::array<T, N>&)>, sizeof...(E)> to_function(const std::tuple<E...>& e);
};

template<class T, size_t N, class... E>
std::array<std::function<T(const std::array<T, N>&)>, sizeof...(E)> to_function(const std::tuple<E...>& e) {
    return _to_array_of_functions::to_array_of_functions<T, N>(e, std::make_index_sequence<sizeof...(E)>{});
}

}

#endif
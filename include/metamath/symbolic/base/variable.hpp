#ifndef SYMBOLIC_VARIABLE_HPP
#define SYMBOLIC_VARIABLE_HPP

#include "integral_constant.hpp"

#include <array>
#include <tuple>
#include <vector>

namespace SYMBOLIC_NAMESPACE {

template<auto N>
struct variable : expression<variable<N>> {
    constexpr operator decltype(N)() const noexcept {
        return N;
    }

    template<class... Args>
    constexpr auto operator()(const std::tuple<Args...>& x) const noexcept {
        static_assert(N < sizeof...(Args));
        return std::get<N>(x);
    }

    template<class... Args>
    constexpr auto operator()(const Args&... args) const noexcept {
        return (*this)(std::forward_as_tuple(args...));
    }

    template<class T, size_t M>
    constexpr T operator()(const std::array<T, M>& x) const noexcept {
        static_assert(N < M);
        return x[N];
    }

    template<class T>
    T operator()(const std::vector<T>& x) const {
        return x[N];
    }

    template<auto X>
    constexpr auto derivative() const noexcept {
        if constexpr (N == X)
            return integral_constant<1>{};
        else
            return integral_constant<0>{};
    }
};

template<uintmax_t N>
constexpr variable<N> simplify(const variable<N>) noexcept {
    return {};
}

}

#endif
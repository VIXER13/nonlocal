#pragma once

#include "integral_constant.hpp"

namespace metamath::symbolic {

template<class T>
class constant : public expression<constant<T>> {
    static_assert(std::is_arithmetic_v<T>);

    const T _value = 0;

public:
    constexpr constant(const T& value) noexcept
        : _value{value} {}

    template<class... Args>
    constexpr const T& operator()(const Args&...) const noexcept {
        return _value;
    }

    template<auto X>
    constexpr integral_constant<0> derivative() const noexcept {
        return {};
    }
};

template<class T>
constexpr constant<T> simplify(const constant<T>& c) noexcept {
    return c;
}

}
#pragma once

#include <concepts>
#include <stdexcept>
#include <string>

namespace metamath::functions {

template<std::integral T, bool Check = true>
constexpr T factorial(const T n) noexcept(!Check) {
    if (Check && n < 0)
        throw std::logic_error{"factorial. Invalid number: " + std::to_string(n)};
    return n > 0 ? n * factorial(n - 1) : 1;
}

}
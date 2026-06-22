#pragma once

#include <metamath/functions/power.hpp>
#include <metamath/utils/operators.hpp>
#include <metamath/utils/constants.hpp>
#include <metamath/types/traits.hpp>

#include <array>
#include <cmath>
#include <numeric>
#include <ranges>

namespace metamath::linear {

template<size_t Exp = 2, std::ranges::random_access_range Container>
auto powered_norm(const Container& x) {
    using T = std::ranges::range_value_t<Container>;
    return std::accumulate(x.begin(), x.end(), T{0}, [](const T sum, const T x) {
        if constexpr (Exp & 1)
            return sum + functions::power<Exp>(std::abs(x));
        return sum + functions::power<Exp>(x);
    });
}

template<std::ranges::random_access_range Container, types::arithmetic Exp>
auto powered_norm(const Container& x, const Exp exp) {
    using T = std::ranges::range_value_t<Container>;
    return std::accumulate(x.begin(), x.end(), T{0}, [exp](const T sum, const T x) {
        return sum + functions::power(std::abs(x), exp);
    });
}

template<size_t Exp = 2, std::ranges::random_access_range Container>
auto norm(const Container& x) {
    using T = std::ranges::range_value_t<Container>;
    if constexpr (Exp == 1)
        return powered_norm<Exp>(x);
    else if constexpr (Exp == 2)
        return std::sqrt(powered_norm<Exp>(x));
    else if constexpr (Exp == 3)
        return std::cbrt(powered_norm<Exp>(x));
    else if constexpr (Exp == constants::Infinity<size_t>)
        return std::abs(*std::max_element(x.begin(), x.end(), 
            [](const T a, const T b) { return std::abs(a) < std::abs(b); }));
    else
        return functions::power(powered_norm<Exp>(x), T{1} / Exp);
}

template<std::ranges::random_access_range Container, types::arithmetic Exp>
auto norm(const Container& x, const Exp exp) {
    using T = std::ranges::range_value_t<Container>;
    return functions::power(powered_norm(x, exp), T{1} / exp);
}

}
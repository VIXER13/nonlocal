#pragma once

#include <metamath/functions/power.hpp>
#include <metamath/utils/constants.hpp>
#include <metamath/types/traits.hpp>

#include <algorithm>
#include <array>
#include <numeric>
#include <ranges>

namespace metamath::linear {

template<size_t Exp = 2, std::ranges::random_access_range Container>
auto powered_norm(const Container& x) {
    using RT = std::ranges::range_value_t<Container>;
    using T = types::container_type_t<RT>;
    static constexpr auto Accumulator = [](const T sum, const T value) {
        return sum + (Exp & 1 ? functions::power<Exp>(std::abs(value)) : functions::power<Exp>(value));
    };
    if constexpr (types::is_array_v<RT>)
        return std::accumulate(x.begin(), x.end(), T{0}, [](const T sum, const RT& x) {
            return std::accumulate(x.begin(), x.end(), sum, Accumulator);
        });
    else
        return std::accumulate(x.begin(), x.end(), T{0}, Accumulator);
}

template<std::ranges::random_access_range Container, types::arithmetic Exp>
auto powered_norm(const Container& x, const Exp exp) {
    using RT = std::ranges::range_value_t<Container>;
    using T = types::container_type_t<RT>;
    const auto accumulator = [exp](const T sum, const T value) {
        return sum + functions::power(std::abs(value), exp);
    };
    if constexpr (types::is_array_v<RT>)
        return std::accumulate(x.begin(), x.end(), T{0}, [&accumulator](const T sum, const RT& x) {
            return std::accumulate(x.begin(), x.end(), sum, accumulator);
        });
    else
        return std::accumulate(x.begin(), x.end(), T{0}, accumulator);
}

template<size_t Exp = 2, std::ranges::random_access_range Container>
auto norm(const Container& x) {
    using RT = std::ranges::range_value_t<Container>;
    using T = types::container_type_t<RT>;
    if constexpr (Exp == 1)
        return powered_norm<Exp>(x);
    else if constexpr (Exp == 2)
        return std::sqrt(powered_norm<Exp>(x));
    else if constexpr (Exp == 3)
        return std::cbrt(powered_norm<Exp>(x));
    else if constexpr (Exp == constants::Infinity<size_t>) {
        T max = T{0};
        static constexpr auto Comprator = [](const T a, const T b) { return std::abs(a) < std::abs(b); };
        if constexpr (types::is_array_v<RT>) {
            for (const auto& value : x)
                if (const T local_max = std::abs(*std::max_element(value.begin(), value.end(), Comprator)); local_max > max)
                    max = local_max;
        } else
            return std::abs(*std::max_element(x.begin(), x.end(), Comprator));
        return max;
    } else
        return functions::power(powered_norm<Exp>(x), T{1} / Exp);
}

template<std::ranges::random_access_range Container, types::arithmetic Exp>
auto norm(const Container& x, const Exp exp) {
    using T = types::container_type_t<std::ranges::range_value_t<Container>>;
    return functions::power(powered_norm(x, exp), T{1} / exp);
}

}
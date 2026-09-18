#pragma once

#include <metamath/types/traits.hpp>

#include <functional>
#include <numeric>
#include <ranges>

namespace metamath::linear {

template<std::ranges::random_access_range Container>
auto scalar_product(const Container& x, const Container& y) {
    using RT = std::ranges::range_value_t<Container>;
    using T = types::container_type_t<RT>;
    if constexpr (!types::is_array_v<Container>)
        if (x.size() != y.size())
            throw std::invalid_argument("Vectors must have the same size for scalar production.");

    if constexpr (types::is_array_v<RT>)
        return std::inner_product(x.begin(), x.end(), y.begin(), T{0}, std::plus<T>(), [](const RT& a, const RT& b) {
            return std::inner_product(a.begin(), a.end(), b.begin(), T{0});
        });
    else
        return std::inner_product(x.begin(), x.end(), y.begin(), T{0});
}

}
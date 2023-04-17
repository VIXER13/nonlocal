#include "init_uniform_ranges.hpp"

#include <exception>
#include <iostream>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> init_uniform_ranges(const size_t size, const size_t count) {
    if (!count)
        throw std::domain_error{"The count parameter cannot be 0!"};
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(count);
    if (const size_t count_per_range = size / count; count_per_range != 0) {
        for(const size_t i : std::ranges::iota_view{0u, count - 1}) {
            const size_t left_bound = count_per_range * i;
            ranges[i] = { left_bound, left_bound + count_per_range };
        }
        ranges.back() = { count_per_range * (count - 1), size };
    } else {
        for(const size_t i : std::ranges::iota_view{0u, size})
            ranges[i] = { i, i + 1 };
        for(const size_t i : std::ranges::iota_view{size, count})
            ranges[i] = { size, size };
    }
    return ranges;
}

}
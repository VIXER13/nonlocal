#include "init_uniform_ranges.hpp"

#include <exception>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> init_uniform_ranges(const size_t size, const size_t count) {
    if (!count)
        throw std::domain_error{"The count parameter cannot be 0!"};
    size_t left_bound = 0u;
    const bool is_small_size = size < count;
    const size_t last = is_small_size ? size : count - 1;
    const size_t count_per_range = is_small_size ? 1u : size / count;
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(count);
    for(const size_t i : std::ranges::iota_view{0u, last}) {
        ranges[i] = { left_bound, left_bound + count_per_range };
        left_bound += count_per_range;
    }
    for(const size_t i : std::ranges::iota_view{last, count})
        ranges[i] = { left_bound, size };
    return ranges;
}

}
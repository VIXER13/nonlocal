#include "uniform_ranges.hpp"

#include <exception>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> uniform_ranges(const size_t size, const size_t ranges_count) {
    if (!ranges_count)
        throw std::domain_error{"The ranges count cannot be 0!"};
    size_t left_bound = 0u;
    const bool is_small_size = size < ranges_count;
    const size_t last = is_small_size ? size : ranges_count - 1;
    const size_t count_per_range = is_small_size ? 1u : size / ranges_count;
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(ranges_count);
    for(const size_t i : std::ranges::iota_view{0u, last}) {
        const size_t right_bound = left_bound + count_per_range;
        ranges[i] = { left_bound, right_bound };
        left_bound = right_bound;
    }
    for(const size_t i : std::ranges::iota_view{last, ranges_count})
        ranges[i] = { left_bound, size };
    return ranges;
}

}
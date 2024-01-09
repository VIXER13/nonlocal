#include "uniform_ranges.hpp"

#include <exception>
#include <numeric>

#include <iostream>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> uniform_ranges(const size_t size, const size_t ranges_count) {
    if (!ranges_count)
        throw std::domain_error{"The ranges count cannot be 0!"};
    size_t left_bound = 0u;
    const bool is_small_size = size < ranges_count;
    const size_t last = is_small_size ? size : ranges_count - 1;
    const size_t per_range = is_small_size ? 1u : size / ranges_count;
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(ranges_count);
    for(const size_t i : std::ranges::iota_view{0u, last}) {
        const size_t right_bound = left_bound + per_range;
        ranges[i] = { left_bound, right_bound };
        left_bound = right_bound;
    }
    for(const size_t i : std::ranges::iota_view{last, ranges_count})
        ranges[i] = { left_bound, size };
    return ranges;
}

std::vector<std::ranges::iota_view<size_t, size_t>> uniform_ranges(const std::vector<size_t> numbers, const size_t ranges_count) {
    if (!ranges_count)
        throw std::domain_error{"The ranges count cannot be 0!"};
    size_t sum = 0u;
    size_t curr_row = 0u;
    size_t curr_range = 0u;
    const size_t mean = std::accumulate(numbers.begin(), numbers.end(), size_t{0}) / ranges_count;
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(ranges_count);
    for(const size_t row : std::ranges::iota_view{0u, numbers.size()}) {
        sum += numbers[row];
        if (sum >= mean) {
            if (curr_range < ranges_count - 1)
                ranges[curr_range] = {curr_row, row};
            else {
                ranges[curr_range] = {curr_range, numbers.size()};
                break;
            }
            curr_row = row;
            ++curr_range;
            sum = 0;
        }
    }
    for(const size_t range : std::ranges::iota_view{curr_range, ranges_count}) {
        ranges[range] = {curr_row, numbers.size()};
        curr_row = numbers.size();
    }
    return ranges;
}

}
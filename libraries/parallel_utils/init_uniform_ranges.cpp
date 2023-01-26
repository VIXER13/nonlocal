#include "init_uniform_ranges.hpp"

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> init_uniform_ranges(const size_t size, const size_t count) {
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(count);
    for(const size_t i : std::ranges::iota_view{0u, count})
        ranges[i] = { size / ranges.size() *  i,
                      size / ranges.size() * (i + 1) + (i == ranges.size() - 1) * size % ranges.size() };
    return ranges;
}

}
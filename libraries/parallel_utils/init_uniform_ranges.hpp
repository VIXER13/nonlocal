#ifndef PARALLEL_UTILS_INIT_UNIFORM_RANGES_HPP
#define PARALLEL_UTILS_INIT_UNIFORM_RANGES_HPP

#include <cstddef>
#include <ranges>
#include <vector>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> init_uniform_ranges(const size_t size, const size_t count);

}

#endif
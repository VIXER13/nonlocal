#ifndef OMP_UTILS_HPP
#define OMP_UTILS_HPP

#include <cstddef>
#include <ranges>
#include <vector>

namespace parallel_utils {

int threads_count();

class OMP_ranges final {
    std::vector<std::ranges::iota_view<size_t, size_t>> _ranges;

public:
    explicit OMP_ranges(const size_t size = 0, const size_t threads = threads_count());

    std::ranges::iota_view<size_t, size_t> get(const size_t thread) const;
    void set(const std::ranges::iota_view<size_t, size_t> range, const size_t thread);
};

}

#endif
#ifndef PARALLEL_UTILS_UNIFORM_RANGES_HPP
#define PARALLEL_UTILS_UNIFORM_RANGES_HPP

#include <Eigen/Sparse>

#include <cstddef>
#include <ranges>
#include <vector>

namespace parallel_utils {

std::vector<std::ranges::iota_view<size_t, size_t>> uniform_ranges(const size_t size, const size_t ranges_count);

// Returns rows ranges where each range has approximately the same elements number.
template<class T, class I>
std::vector<std::ranges::iota_view<size_t, size_t>> uniform_ranges(const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& matrix, const size_t ranges_count) {
    if (!ranges_count)
        throw std::domain_error{"The ranges count cannot be 0!"};
    size_t curr_row = 0u;
    size_t curr_range = 0u;
    const size_t mean = matrix.nonZeros() / ranges_count;
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(ranges_count);
    for(const size_t row : std::ranges::iota_view{0u, size_t(matrix.rows())})
        if (const size_t sum = matrix.outerIndexPtr()[row + 1] - matrix.outerIndexPtr()[curr_row]; sum > mean) {
            if (curr_range < ranges_count - 1)
                ranges[curr_range] = {curr_row, row};
            else {
                ranges[curr_range] = {curr_row, size_t(matrix.rows())};
                break;
            }
            curr_row = row;
            ++curr_range;
        }
    for (const size_t range : std::ranges::iota_view{curr_range, ranges_count}) {
        ranges[range] = {curr_row, size_t(matrix.rows())};
        curr_row = size_t(matrix.rows());
    }
    return ranges;
}

}

#endif
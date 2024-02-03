#include "MPI_utils.hpp"

#include "uniform_ranges.hpp"

#include <iostream>

namespace parallel {

int MPI_rank() {
    int rank = 0;
#if MPI_BUILD
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
}

int MPI_size() {
    int size = 1;
#if MPI_BUILD
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    return size;
}

bool is_last_process() {
    return MPI_rank() == MPI_size() - 1;
}

std::vector<std::ranges::iota_view<size_t, size_t>> rows_distribution(const size_t rows) {
    std::vector<std::ranges::iota_view<size_t, size_t>> ranges(MPI_size());
    std::ranges::iota_view<size_t, size_t> range = {0, rows};
#if MPI_BUILD
    MPI_Allgather(&range, sizeof(range), MPI_BYTE, ranges.data(), sizeof(range), MPI_BYTE, MPI_COMM_WORLD);
    for(const size_t i : std::ranges::iota_view{1u, ranges.size()})
        ranges[i] = {*ranges[i - 1].end(), *ranges[i - 1].end() + ranges[i].size()};
#else
    ranges = {range};
#endif
    return ranges;
}

MPI_ranges::MPI_ranges(const size_t size)
    : _ranges{uniform_ranges(size, MPI_size())} {}

MPI_ranges::MPI_ranges(const std::vector<std::ranges::iota_view<size_t, size_t>>& ranges)
    : _ranges{ranges} {
    if (_ranges.size() != MPI_size())
        throw std::domain_error{"The ranges count when initializing MPI_ranges must match the running MPI processes count."};
}

std::ranges::iota_view<size_t, size_t> MPI_ranges::get(const size_t process) const {
    return _ranges[process];
}

void MPI_ranges::set(const std::ranges::iota_view<size_t, size_t> range, const size_t process) {
    _ranges[process] = range;
}

}
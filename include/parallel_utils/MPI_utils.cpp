#include "MPI_utils.hpp"

namespace parallel_utils {

int MPI_rank() noexcept {
    int rank = 0;
#if MPI_USE
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
}

int MPI_size() noexcept {
    int size = 1;
#if MPI_USE
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    return size;
}

MPI_ranges::MPI_ranges()
    : _ranges(MPI_size(), std::array<size_t, 2>{}) {}

const std::vector<std::array<size_t, 2>>& MPI_ranges::ranges() const noexcept {
    return _ranges;
}

const std::array<size_t, 2>& MPI_ranges::range(const size_t rank) const noexcept {
    return _ranges[rank];
}

std::array<size_t, 2>& MPI_ranges::range(const size_t rank) noexcept {
    return _ranges[rank];
}

}
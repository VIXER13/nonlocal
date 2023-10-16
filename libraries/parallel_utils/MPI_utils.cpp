#include "MPI_utils.hpp"

#include "init_uniform_ranges.hpp"

namespace parallel_utils {

int MPI_rank() {
    int rank = 0;
#if MPI_USED
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
}

int MPI_size() {
    int size = 1;
#if MPI_USED
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    return size;
}

bool is_last_process() {
    return MPI_rank() == MPI_size() - 1;
}

MPI_ranges::MPI_ranges(const size_t size)
    : _ranges{init_uniform_ranges(size, MPI_size())} {}

std::ranges::iota_view<size_t, size_t> MPI_ranges::get(const size_t process) const {
    return _ranges[process];
}

void MPI_ranges::set(const std::ranges::iota_view<size_t, size_t> range, const size_t process) {
    _ranges[process] = range;
}

}
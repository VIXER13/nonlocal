#include "MPI_utils.hpp"

namespace parallel_utils {

int MPI_rank() {
    int rank = 0;
#if MPI_USE
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    return rank;
}

int MPI_size() {
    int size = 1;
#if MPI_USE
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
    return size;
}

MPI_ranges::MPI_ranges(const size_t size)
    : _ranges(MPI_size()) {
    for(const size_t process : std::ranges::iota_view{0u, _ranges.size()})
        _ranges[process] = { size / _ranges.size() *  process,
                             size / _ranges.size() * (process + 1) + (process == _ranges.size() - 1) * size % _ranges.size() };
}

std::ranges::iota_view<size_t, size_t> MPI_ranges::get(const size_t process) const {
    return _ranges[rank];
}

void MPI_ranges::set(const std::ranges::iota_view<size_t, size_t> range, const size_t process) {
    _ranges[rank] = range;
}

}
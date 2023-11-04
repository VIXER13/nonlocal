#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#if MPI_BUILD
#   include <mpi.h>
#endif

#include <array>
#include <vector>
#include <ranges>

#include <stddef.h>

namespace parallel_utils {

int MPI_rank();
int MPI_size();
bool is_last_process();
std::vector<std::ranges::iota_view<size_t, size_t>> rows_distribution(const size_t rows);

class MPI_ranges final {
    std::vector<std::ranges::iota_view<size_t, size_t>> _ranges;

public:
    explicit MPI_ranges(const size_t size = 0);
    explicit MPI_ranges(const std::vector<std::ranges::iota_view<size_t, size_t>>& ranges);

    std::ranges::iota_view<size_t, size_t> get(const size_t process = MPI_rank()) const;
    void set(const std::ranges::iota_view<size_t, size_t> range, const size_t process = MPI_rank());
};

template<class U, class Vector>
Vector all_to_all(const Vector& sendbuf, const MPI_ranges& ranges) {
// #if MPI_BUILD
//     std::vector<int> sendcounts(MPI_size(),        sizeof(U) * (ranges.range().back() - ranges.range().front())),
//                      sdispls   (sendcounts.size(), sizeof(U) *  ranges.range().front()),
//                      recvcounts(sendcounts.size()), rdispls(sendcounts.size());
//     for(const size_t i : std::ranges::iota_view{size_t{0}, sendcounts.size()}) {
//         recvcounts[i] = sizeof(U) * (ranges.range(i).back() - ranges.range(i).front());
//         rdispls[i]    = sizeof(U) *  ranges.range(i).front();
//     }
//     Vector recvbuf(sendbuf.size()); // for old version mpich, when strict sendbuf and recvbuf don't support
//     // double cast also for old mpich version
//     MPI_Alltoallv(const_cast<void*>(static_cast<const void*>(sendbuf.data())), sendcounts.data(), sdispls.data(), MPI_BYTE,
//                   recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_BYTE, MPI_COMM_WORLD);
// #else
    Vector recvbuf = sendbuf;
// #endif
    return recvbuf;
}

template<class T>
constexpr std::enable_if_t<std::is_floating_point_v<T>, T> reduce(T local_sum) {
// #if MPI_BUILD
//     T global_sum = T{0};
//     MPI_Reduce(&local_sum, &global_sum, 1, std::is_same_v<T, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//     return global_sum;
// #else
    return local_sum;
// #endif
}

}

#endif
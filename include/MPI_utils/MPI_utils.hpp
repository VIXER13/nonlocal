#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <array>
#include <vector>
#if MPI_USE
    #include <mpi.h>
#endif

namespace MPI_utils {

int MPI_rank() noexcept;
int MPI_size() noexcept;

class MPI_ranges final {
    std::vector<std::array<size_t, 2>> _ranges;

public:
    MPI_ranges();

    const std::vector<std::array<size_t, 2>>& ranges() const noexcept;
    const std::array<size_t, 2>& range(const size_t rank = MPI_rank()) const noexcept;
    std::array<size_t, 2>& range(const size_t rank = MPI_rank()) noexcept;
};

template<class U, class Vector>
Vector all_to_all(const Vector& sendbuf, const MPI_ranges& ranges) {
    std::vector<int> sendcounts(MPI_size(),        sizeof(U) * (ranges.range().back() - ranges.range().front())),
                     sdispls   (sendcounts.size(), sizeof(U) *  ranges.range().front()),
                     recvcounts(sendcounts.size()), rdispls(sendcounts.size());
    for(size_t i = 0; i < recvcounts.size(); ++i) {
        recvcounts[i] = sizeof(U) * (ranges.range(i).back() - ranges.range(i).front());
        rdispls[i]    = sizeof(U) *  ranges.range(i).front();
    }
#if MPI_USE
    Vector recvbuf(sendbuf.size()); // for old version mpich, when strict sendbuf and recvbuf don't support
    // double cast also for old mpich version
    MPI_Alltoallv(const_cast<void*>(static_cast<const void*>(sendbuf.data())), sendcounts.data(), sdispls.data(), MPI_BYTE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_BYTE, MPI_COMM_WORLD);
#else
    Vector recvbuf = sendbuf;
#endif
    return recvbuf;
}

template<class T>
constexpr std::enable_if_t<std::is_floating_point_v<T>, T> reduce(T local_sum) {
#if MPI_USE
    T global_sum = T{0};
    MPI_Reduce(&local_sum, &global_sum, 1, std::is_same_v<T, float> ? MPI_FLOAT : MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum;
#else
    return local_sum;
#endif
}

}

#endif
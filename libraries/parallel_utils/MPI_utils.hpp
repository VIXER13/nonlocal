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

template<class Vector>
Vector all_to_all(const Vector& sendbuf, const MPI_ranges& ranges) {
#if MPI_BUILD
    static constexpr size_t data_size = sizeof(sendbuf[0]);
    std::vector<int> sendcounts(MPI_size(),        data_size * ranges.get().size()),
                     sdispls   (sendcounts.size(), data_size * ranges.get().front()),
                     recvcounts(sendcounts.size()), rdispls(sendcounts.size());
    for(const size_t i : std::ranges::iota_view{0u, sendcounts.size()}) {
        recvcounts[i] = data_size * ranges.get(i).size();
        rdispls[i]    = data_size * ranges.get(i).front();
    }
    Vector recvbuf(sendbuf.size()); // for old version mpich, when strict sendbuf and recvbuf don't support
    // double cast also for old mpich version
    MPI_Alltoallv(const_cast<void*>(static_cast<const void*>(sendbuf.data())), sendcounts.data(), sdispls.data(), MPI_BYTE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_BYTE, MPI_COMM_WORLD);
#else
    const Vector recvbuf = sendbuf;
#endif
    return recvbuf;
}

template<class Recv, class Send>
void reduce_vector(Recv& recvbuf, const Send& sendbuf) {
#if MPI_BUILD
    if (recvbuf.size() != sendbuf.size())
        throw std::logic_error{"Cannot perform reduction, recvbuf and sendbuf are of different sizes."};
    if (recvbuf.data() == sendbuf.data())
        throw std::domain_error{"During reduction, the send and receive buffers must not overlap."};
    static constexpr bool is_float = std::is_same_v<std::remove_cvref_t<decltype(sendbuf[0])>, float>;
    MPI_Allreduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), is_float ? MPI_FLOAT : MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    recvbuf = sendbuf;
#endif
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
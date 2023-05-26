#ifndef METAMATH_UNIFORM_PARTITION_HPP
#define METAMATH_UNIFORM_PARTITION_HPP

#include <array>

namespace metamath::utils {

class _uniform_partition final {
    constexpr explicit _uniform_partition() noexcept = default;

    template<size_t N, class T, size_t... I>
    static constexpr std::array<T, N> uniform_partition(const std::array<T, 2>& segment, const std::index_sequence<I...>) noexcept {
        static_assert(N > 0, "Unable to create 0 nodes on a segment.");
        if constexpr (N == 1)
            return {(segment.front() + segment.back()) / T(2)};
        const T step = (segment.back() - segment.front()) / T(N - 1);
        return {segment.front() + step * T(I)...};
    }

public:
    template<size_t N, class T>
    friend constexpr std::array<T, N> uniform_partition(const std::array<T, 2>& segment) noexcept;
};

template<size_t N, class T>
constexpr std::array<T, N> uniform_partition(const std::array<T, 2>& segment) noexcept {
    return _uniform_partition::uniform_partition<N>(segment, std::make_index_sequence<N>{});
}

}

#endif
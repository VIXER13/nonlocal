#ifndef MESH_UTILS_HPP
#define MESH_UTILS_HPP

#include "metamath.hpp"
#include <array>

namespace nonlocal::mesh::utils {

template<class T, size_t N>
T distance(const std::array<T, N>& A, const std::array<T, N>& B) noexcept {
    T sum = 0;
    for(size_t i = 0; i < N; ++i)
        sum += metamath::function::power<2>(A[i] - B[i]);
    return std::sqrt(sum);
}

}

#endif
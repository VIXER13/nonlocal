#pragma once

#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>

#include <concepts>
#include <functional>
#include <vector>

namespace nonlocal::unit_tests {

template<std::floating_point T>
T max_error(const std::vector<T>& x, const std::vector<T>& y) {
    using namespace metamath::operators;
    static constexpr auto Inf = metamath::constants::Infinity<size_t>;
    return metamath::linear::norm<Inf>(x - y);
}

template<class Vector, std::floating_point T, std::integral I, class Expected>
T norm_error(const Vector& actual, const mesh::mesh_container_2d<T, I>& mesh, const Expected& function) {
    using metamath::linear::norm;
    using namespace metamath::operators;
    static constexpr auto Inf = metamath::constants::Infinity<size_t>;
    const auto discrete_function = nonlocal::mesh::utils::discrete(mesh, function);
    return norm<Inf>(actual - discrete_function) / norm<Inf>(discrete_function);
}

template<std::floating_point T>
T L2_norm(const std::vector<T>& x, const std::vector<T>& y) {
    if (x.size() != y.size()) throw std::invalid_argument("Vectors must be the same size.");
    T norm = std::transform_reduce(x.begin(), x.end(), y.begin(), static_cast<T>(0.0), std::plus<>(),
    [](T x, T y) {
        T diff = x - y;
        return diff * diff;
    });
    return std::sqrt(norm);
}


}
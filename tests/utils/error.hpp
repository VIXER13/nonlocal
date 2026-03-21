#pragma once

#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>

#include <concepts>
#include <functional>
#include <vector>

namespace nonlocal::unit_tests {

template<std::floating_point T>
T max_norm(const std::vector<T>& x) {
    return std::abs(*std::max_element(x.begin(), x.end(), [](const T a, const T b) { return std::abs(a) < std::abs(b); }));
}

template<std::floating_point T>
T max_error(const std::vector<T>& x, const std::vector<T>& y) {
    using namespace metamath::functions;
    return max_norm(x - y);
}

template<std::floating_point T, std::integral I, class Expected>
T norm_error(const std::vector<T>& actual, const mesh::mesh_container_2d<T, I>& mesh, const Expected& function) {
    const auto discrete_function = nonlocal::mesh::utils::discrete(mesh, function);
    return max_error(actual, discrete_function) / max_norm(discrete_function);
}

}
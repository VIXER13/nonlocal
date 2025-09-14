#pragma once

#include <mesh/mesh_1d/mesh_1d.hpp>

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

}
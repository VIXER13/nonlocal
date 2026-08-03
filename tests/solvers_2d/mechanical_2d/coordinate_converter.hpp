#pragma once

#include <array>
#include <concepts>
#include <cmath>

namespace nonlocal::unit_tests {

enum : size_t {RR, FF, RF};

template<std::floating_point T>
std::array<T, 2> polar_to_cartesian(const T value, const std::array<T, 2>& point) noexcept {
    const T angle = std::atan2(point[1], point[0]);
    return {value * std::cos(angle), value * std::sin(angle)};
}

template<std::floating_point T>
std::array<T, 3> polar_to_cartesian(const std::array<T, 3>& polar_tensor, const std::array<T, 2>& point) noexcept {
    const T angle = std::atan2(point[1], point[0]);
    const T cos = std::cos(angle);
    const T sin = std::sin(angle);
    const T cos2 = cos * cos;
    const T sin2 = sin * sin;
    const T shear = 2 * polar_tensor[RF] * sin * cos;
    return {
        polar_tensor[RR] * cos2 + polar_tensor[FF] * sin2 - shear,
        polar_tensor[RR] * sin2 + polar_tensor[FF] * cos2 + shear,
        (polar_tensor[RR] - polar_tensor[FF]) * sin * cos + polar_tensor[RF] * (cos2 - sin2)
    };
}

}
#pragma once

#include <cstdint>
#include <concepts>

namespace nonlocal {

enum axis : uint8_t { X, Y, Z };
enum : uint8_t { XX, YY, XY, YX = XY };

enum class physics_t : uint8_t {
    THERMAL,
    MECHANICAL
};

enum class theory_t : bool {
    LOCAL,
    NONLOCAL
};

template<std::floating_point T>
constexpr T Nonlocal_Threshold = T{0.999};

template<std::floating_point T>
constexpr theory_t theory_type(const T local_weight) noexcept {
    return local_weight < Nonlocal_Threshold<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
}

template<std::floating_point T>
constexpr bool is_local(const T local_weight) noexcept {
    return theory_type(local_weight) == theory_t::LOCAL;
}

template<std::floating_point T>
constexpr bool is_nonlocal(const T local_weight) noexcept {
    return theory_type(local_weight) == theory_t::NONLOCAL;
}

template<std::floating_point T>
constexpr T nonlocal_weight(const T local_weight) noexcept {
    return T{1} - local_weight;
}

template<std::floating_point T>
inline constexpr T NEUMANN_PROBLEM_ERROR_THRESHOLD = std::is_same_v<T, float> ? T{1e-5} : T{1e-10};

}
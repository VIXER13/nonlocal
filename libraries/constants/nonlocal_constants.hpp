#pragma once

#include <cstdint>
#include <concepts>

namespace nonlocal {

inline constexpr auto EMPTY_FUNCTION = []() constexpr noexcept {};

enum axis : uint8_t {
    X,
    Y,
    Z
};

enum class material_t : uint8_t {
    ISOTROPIC,
    ORTHOTROPIC,
    ANISOTROPIC
};

enum class physics_t : uint8_t {
    THERMAL,
    MECHANICAL
};

enum class coefficients_t : uint8_t {
    CONSTANTS,
    SPACE_DEPENDENT,
    SOLUTION_DEPENDENT
};

enum class boundary_condition_t : uint8_t {
    FIRST_KIND = 1,
    SECOND_KIND,
    THIRD_KIND,
    FOURTH_KIND,
    FIFTH_KIND
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
constexpr T nonlocal_weight(const T local_weight) noexcept {
    return T{1} - local_weight;
}

template<std::floating_point T>
inline constexpr T Stefan_Boltzmann_Constant = T{5.67036713e-8};

template<std::floating_point T>
inline constexpr T NEUMANN_PROBLEM_ERROR_THRESHOLD = std::is_same_v<T, float> ? T{1e-5} : T{1e-10};

}
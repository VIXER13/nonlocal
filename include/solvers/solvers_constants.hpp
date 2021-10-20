#ifndef BASE_CONSTANTS_HPP
#define BASE_CONSTANTS_HPP

#include <cstdint>

namespace nonlocal {

enum class boundary_condition_t : uint8_t {
    FIRST_KIND,
    SECOND_KIND,
    THIRD_KIND,
    FOURTH_KIND
};

enum class material_t : uint8_t {
    ISOTROPIC,
    ORTHOTROPIC,
    ANISOTROPIC
};

enum class theory_t : uint8_t {
    LOCAL,
    NONLOCAL
};

template<class T>
inline constexpr T MAX_NONLOCAL_WEIGHT = T{0.999};

template<class T>
inline constexpr T MAX_NEUMANN_PROBLEM_ERROR = T{1e-5};

}

#endif
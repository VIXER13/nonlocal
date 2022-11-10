#ifndef NONLOCAL_CONSTANTS_HPP
#define NONLOCAL_CONSTANTS_HPP

#include <cstdint>
#include <concepts>

namespace nonlocal {

enum axis : uint8_t {
    X,
    Y,
    Z
};

enum class theory_t : bool {
    LOCAL,
    NONLOCAL
};

template<std::floating_point T>
inline constexpr T MAX_NONLOCAL_WEIGHT = T{0.999};

template<std::floating_point T>
constexpr theory_t theory_type(const T local_weight) noexcept {
    return local_weight < MAX_NONLOCAL_WEIGHT<T> ? theory_t::NONLOCAL : theory_t::LOCAL;
}

template<std::floating_point T>
constexpr T nonlocal_weight(const T local_weight) noexcept {
    return T{1} - local_weight;
}

}

#endif
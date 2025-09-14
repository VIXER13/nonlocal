#pragma once

#include <concepts>
#include <limits>

namespace metamath::constants {

template<std::integral T>
inline constexpr T Infinity = std::numeric_limits<T>::max();

template<std::floating_point T>
inline constexpr T Stefan_Boltzmann_Constant = T{5.67036713e-8};

}
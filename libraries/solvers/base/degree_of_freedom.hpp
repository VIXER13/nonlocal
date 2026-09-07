#pragma once

#include <array>

namespace nonlocal {

template<class T>
struct degree_of_freedom : std::integral_constant<size_t, 1> {};

template<class T, size_t Dimension>
struct degree_of_freedom<std::array<T, Dimension>> : std::integral_constant<size_t, Dimension> {};

template<class T>
constexpr size_t DoF = degree_of_freedom<T>::value;

}
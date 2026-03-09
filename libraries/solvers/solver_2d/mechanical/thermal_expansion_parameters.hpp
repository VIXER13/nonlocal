#pragma once

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::mechanical {

// single constant
template<class T>
using isotropic_thermal_expansion_t = T;
// a[0]  0
//  0   a[1]
template<class T>
using orthotropic_thermal_expansion_t = std::array<T, 2>;
// a[0] a[2]
// a[2] a[1]
template<class T>
using anisotropic_thermal_expansion_t = std::array<T, 3>;
template<class T>
using thermal_expansion_t = std::variant<
    isotropic_thermal_expansion_t<T>,
    orthotropic_thermal_expansion_t<T>,
    anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
using raw_isotropic_thermal_expansion_t = isotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_orthotropic_thermal_expansion_t = orthotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_anisotropic_thermal_expansion_t = anisotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_thermal_expansion_t = std::variant<
    raw_isotropic_thermal_expansion_t<T>,
    raw_orthotropic_thermal_expansion_t<T>,
    raw_anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
using evaluated_isotropic_thermal_expansion_t = evaluated_parameters<isotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_orthotropic_thermal_expansion_t = evaluated_parameters<orthotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_anisotropic_thermal_expansion_t = evaluated_parameters<anisotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_thermal_expansion_t = std::variant<
    evaluated_isotropic_thermal_expansion_t<T>,
    evaluated_orthotropic_thermal_expansion_t<T>,
    evaluated_anisotropic_thermal_expansion_t<T>
>;

}
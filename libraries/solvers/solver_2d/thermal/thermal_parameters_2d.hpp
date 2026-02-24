#pragma once

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::thermal {

// single constant
template<class T>
using isotropic_conductivity_t = T;
// l[0]  0
//  0   l[1]
template<class T>
using orthotropic_conductivity_t = std::array<T, 2>;
// l[0] l[2]
// l[2] l[1]
template<class T>
using anisotropic_conductivity_t = std::array<T, 3>;
template<class T>
using conductivity_t = std::variant<
    isotropic_conductivity_t<T>,
    orthotropic_conductivity_t<T>,
    anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
using raw_isotropic_conductivity_t = isotropic_conductivity_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_orthotropic_conductivity_t = orthotropic_conductivity_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_anisotropic_conductivity_t = anisotropic_conductivity_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_conductivity_t = std::variant<
    raw_isotropic_conductivity_t<T>,
    raw_orthotropic_conductivity_t<T>,
    raw_anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
using evaluated_isotropic_conductivity_t = evaluated_parameters<isotropic_conductivity_t<T>>;
template<std::floating_point T>
using evaluated_orthotropic_conductivity_t = evaluated_parameters<orthotropic_conductivity_t<T>>;
template<std::floating_point T>
using evaluated_anisotropic_conductivity_t = evaluated_parameters<anisotropic_conductivity_t<T>>;
template<std::floating_point T>
using evaluated_conductivity_t = std::variant<
    evaluated_isotropic_conductivity_t<T>,
    evaluated_orthotropic_conductivity_t<T>,
    evaluated_anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
struct parameter_2d final {
    raw_conductivity_t<T> conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
    T relaxation_time = T{0};
};

template<std::floating_point T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

template<std::floating_point T>
using evaluated_conductivity_2d = std::unordered_map<std::string, equation_parameters<2, T, evaluated_conductivity_t>>;

}
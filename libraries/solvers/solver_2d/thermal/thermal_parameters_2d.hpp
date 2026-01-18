#pragma once

#include <metamath/metamath.hpp>
#include <solvers/base/equation_parameters.hpp>

#include <string>
#include <unordered_map>
#include <variant>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
using isotropic_conductivity_t = coefficient_t<T, 2>;

template<std::floating_point T>
using orthotropic_conductivity_t = std::array<coefficient_t<T, 2>, 2>;

template<std::floating_point T>
using anisotropic_conductivity_t = std::array<coefficient_t<T, 2>, 3>;

template<std::floating_point T>
using conductivity_t = std::variant<
    isotropic_conductivity_t<T>,
    orthotropic_conductivity_t<T>,
    anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
struct parameter_2d final {
    conductivity_t<T> conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
    T relaxation_time = T{0};
};

template<std::floating_point T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

template<std::floating_point T>
using evaluated_isotropic_conductivity_t = std::variant<T, metamath::types::vector_with_shifted_index<T>>;

template<std::floating_point T>
using evaluated_orthotropic_conductivity_t = std::variant<
    std::array<T, 2>,
    metamath::types::vector_with_shifted_index<std::array<T, 2>>
>;

template<std::floating_point T>
using evaluated_anisotropic_conductivity_t = std::variant<
    std::array<T, 3>,
    metamath::types::vector_with_shifted_index<std::array<T, 3>>
>;

template<std::floating_point T>
using evaluated_conductivity_t = std::variant<
    evaluated_isotropic_conductivity_t<T>,
    evaluated_orthotropic_conductivity_t<T>,
    evaluated_anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
using evaluated_conductivity_2d = std::unordered_map<std::string, equation_parameters<2, T, evaluated_conductivity_t>>;

}
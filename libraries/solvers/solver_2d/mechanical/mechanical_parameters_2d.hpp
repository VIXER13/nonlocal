#pragma once

#include "elastic_parameters.hpp"
#include "thermal_expansion_parameters.hpp"

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
struct raw_mechanical_parameters_t final {
    elastic_parameters_t<T> elastic;
    raw_thermal_expansion_t<T> thermal_expansion;
};

template<std::floating_point T>
struct evaluated_mechanical_parameters_t final {
    evaluated_hook_matrix_t<T> elastic;
    evaluated_thermal_expansion_t<T> thermal_expansion;
};

template<std::floating_point T>
using raw_mechanical_parameters = std::unordered_map<std::string, equation_parameters<2, T, raw_mechanical_parameters_t>>;

template<std::floating_point T>
using evaluated_mechanical_parameters = std::unordered_map<std::string, equation_parameters<2, T, evaluated_mechanical_parameters_t>>;

}
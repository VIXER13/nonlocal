#pragma once

#include "mechanical_boundary_conditions_1d.hpp"
#include "mechanical_parameters_1d.hpp"

#include <solvers/solver_1d/base/problem_settings.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
problem_settings init_problem_settings(const parameters_1d<T>& parameters,
                                       const mechanical_boundaries_conditions_1d<T>& boundaries_conditions,
                                       const bool is_stationary_problem) {
    static constexpr auto is_force = [](const auto& condition) noexcept {
        return  bool(dynamic_cast<const pressure_1d<T>*>(condition.get())) &&
               !bool(dynamic_cast<const combined_loading_1d<T>*>(condition.get()));
    };
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) noexcept {
        return !is_constant<T, 1>(parameter.physical.youngs_modulus) ||
               !is_constant<T, 1>(parameter.physical.density);
    };
    static constexpr auto is_solution_dependent_parameters = [](const auto& parameter) noexcept {
        return std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.youngs_modulus) ||
               std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.density);
    };
    return {
        .theories = theories_types(parameters),
        .is_neumann = is_stationary_problem ? std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_force) : false,
        .is_nonlinear_boundary = false,
        .is_nonconstant_parameters = std::any_of(parameters.begin(), parameters.end(), is_nonconstant_parameters),
        .is_solution_dependent = std::any_of(parameters.begin(), parameters.end(), is_solution_dependent_parameters),
        .is_first_kind = { 
            bool(dynamic_cast<displacement_1d<T>*>(boundaries_conditions.front().get())),
            bool(dynamic_cast<displacement_1d<T>*>(boundaries_conditions.back ().get()))
        }
    };
}

}
#pragma once

#include "mechanical_boundary_conditions_1d.hpp"
#include "mechanical_parameters_1d.hpp"

#include <solvers/solver_1d/base/problem_settings.hpp>

namespace nonlocal::solver_1d::mechanical {

template<class T>
problem_settings init_problem_settings(const parameters_1d<T>& parameters,
                                       const boundaries_data_1d<T>& boundaries,
                                       const bool is_stationary_problem) {
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) noexcept {
        return !is_constant<T, 1>(parameter.physical.youngs_modulus) || !is_constant<T, 1>(parameter.physical.density);
    };
    static constexpr auto is_solution_dependent = [](const auto& parameter) noexcept {
        return std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.youngs_modulus) ||
               std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.density);
    };
    return init_boundary_settings(
        problem_settings{
            .theories = theories_types(parameters),
            .is_nonconstant_parameters = std::any_of(parameters.begin(), parameters.end(), is_nonconstant_parameters),
            .is_solution_dependent = std::any_of(parameters.begin(), parameters.end(), is_solution_dependent)
        },
        boundaries, is_stationary_problem);
}

template<class T>
problem_settings init_problem_settings(const parameters_1d<T>& parameters,
                                       const mechanical_boundaries_conditions_1d<T>& boundaries,
                                       const bool is_stationary_problem) {
    return init_problem_settings(parameters, normalize_boundaries_1d(boundaries), is_stationary_problem);
}

} // namespace nonlocal::solver_1d::mechanical

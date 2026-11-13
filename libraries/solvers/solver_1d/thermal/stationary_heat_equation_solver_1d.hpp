#pragma once

#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "init_problem_settings.hpp"
#include "radiation_condition_1d.hpp"

#include <solvers/solver_1d/base/stationary_solver_1d.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
struct stationary_equation_parameters_1d final {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T energy = T{ 0 };
};

template<class T>
stationary_solver_parameters_1d<T> to_solver_parameters(const stationary_equation_parameters_1d<T>& parameters) {
    return { .right_part = parameters.right_part,
             .initial_distribution = parameters.initial_distribution,
             .tolerance = parameters.tolerance,
             .max_iterations = parameters.max_iterations,
             .integral_value = parameters.energy };
}

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_right_part(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                    const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                    const stationary_equation_parameters_1d<T>& additional_parameters,
                                                    const bool is_neumann = false) {
    return init_stationary_right_part_1d(*mesh, boundaries_conditions, to_solver_parameters(additional_parameters), is_neumann);
}

template<class T, class I>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const parameters_1d<T>& parameters,
    const thermal_boundaries_conditions_1d<T>& boundaries,
    const stationary_equation_parameters_1d<T>& additional_parameters) {
    const auto settings = init_problem_settings(parameters, boundaries, true);
    const auto solution = stationary_solver_1d<T, I, physics_t::THERMAL>(
        mesh, parameters, boundaries, settings, to_solver_parameters(additional_parameters), &parameter_1d<T>::conductivity,
        convection_condition_1d<T, I>, radiation_condition_1d<true, T, I>);
    return heat_equation_solution_1d<T>{ mesh, parameters, solution };
}

} // namespace nonlocal::solver_1d::thermal

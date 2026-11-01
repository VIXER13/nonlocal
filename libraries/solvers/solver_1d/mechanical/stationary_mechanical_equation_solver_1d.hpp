#pragma once

#include "init_problem_settings.hpp"
#include "mechanical_equation_solution_1d.hpp"

#include <solvers/solver_1d/base/stationary_solver_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
using stationary_equation_parameters_1d = stationary_solver_parameters_1d<T>;

template<std::floating_point T>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_right_part(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                    const boundaries_data_1d<T>& boundaries_conditions,
                                                    const stationary_equation_parameters_1d<T>& additional_parameters,
                                                    const bool is_neumann = false) {
    return init_stationary_right_part_1d(*mesh, boundaries_conditions, additional_parameters, is_neumann);
}

template<std::floating_point T, std::integral I>
mechanical_equation_solution_1d<T> stationary_mechanical_equation_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const parameters_1d<T>& parameters,
    const mechanical_boundaries_conditions_1d<T>& physical_boundaries,
    const stationary_equation_parameters_1d<T>& additional_parameters) {
    const auto boundaries = normalize_boundaries_1d(physical_boundaries);
    const auto settings = init_problem_settings(parameters, boundaries, true);
    const auto solution = stationary_solver_1d<T, I>(
        mesh, parameters, boundaries, settings, additional_parameters, &parameter_1d<T>::youngs_modulus);
    return mechanical_equation_solution_1d<T>{ mesh, parameters, solution };
}

} // namespace nonlocal::solver_1d::mechanical

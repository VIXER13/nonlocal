#pragma once

#include "evaluate_mechanical_parameters.hpp"
#include "stiffness_matrix_2d.hpp"
#include "mechanical_boundary_conditions_2d.hpp"
#include "mechanical_solution_2d.hpp"
#include "init_problem_settings.hpp"
#include "temperature_condition_2d.hpp"

#include <solvers/base/utils.hpp>
#include <solvers/slae/init_solver.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

#include <optional>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
mechanical::mechanical_solution_2d<T> equilibrium_equation(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                           const raw_mechanical_parameters<T>& parameters,
                                                           const mechanical_boundaries_conditions_2d<T>& boundaries_conditions,
                                                           const std::vector<T>& delta_temperature = {},
                                                           const std::function<std::array<T, 2>(const std::array<T, 2>&)>& right_part = nullptr) {
    const auto settings = init_problem_settings(mesh->container(), parameters, boundaries_conditions);
    log_problem_settings(settings);
    const auto evaluated_parameters = evaluate_mechanical_parameters(*mesh, parameters, delta_temperature);

    stiffness_matrix<T> stiffness{*mesh};
    stiffness.compute(evaluated_parameters, settings);
    std::vector<std::array<T, 2>> f(mesh->container().nodes_count(), std::array<T, 2>{});
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    if (right_part)
        integrate_right_part(f, *mesh, right_part);
    // temperature_condition(f, *mesh, evaluated_parameters);
    // boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, stiffness.matrix().bound());

    auto solver = slae::init_iterative_solver(stiffness.matrix(), settings.is_symmetric());
    // if (settings.is_nonlocal()) {
    //     const auto theroires_setter = std::views::all(mesh.container().groups_2d()) |
    //                               std::views::transform([](const std::string& group) { return std::pair{group, theory_t::LOCAL}; });
    //     settings.theories = std::unordered_map<std::string, theory_t>(theroires_setter.begin(), theroires_setter.end());
    //     stiffness_matrix<T> local_stiffness{mesh};
    //     local_stiffness.nodes_for_processing(std::ranges::iota_view{0zu, mesh->container().nodes_count()});
    //     local_stiffness.compute(evaluated_parameters, settings);
    //     if (auto preconditioner = slae::init_preconditioner(std::move(local_stiffness.matrix()), settings.is_symmetric()))
    //         solver->preconditioner(std::move(preconditioner));
    //     else
    //         logger::warning() << "The preconditioner could not be calculated, "
    //                           << "the preconditioner was switched to the Identity." << std::endl;
    // }
    auto solution = mechanical_solution_2d{mesh, evaluated_parameters, solver->solve(f)};
    solution.calc_strain_and_stress(evaluated_parameters);
    return solution;
}

}
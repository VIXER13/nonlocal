#pragma once

#include "conductivity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "radiation_condition_2d.hpp"
#include "init_problem_settings.hpp"
#include "heat_equation_solution_2d.hpp"

#include <solvers/base/utils.hpp>
#include <solvers/slae/init_solver.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T>
struct stationary_equation_parameters_2d final {
    std::function<T(const std::array<T, 2>&)> right_part;
    std::function<T(const std::array<T, 2>&)> initial_distribution;
    T energy = T{0}; // integral condition for Neumann problem
    T tolerance = 100 * std::numeric_limits<T>::epsilon();
    size_t max_iterations = 100;
};

template<std::floating_point T>
std::vector<T> init_right_part(const mesh::mesh_2d<T>& mesh,
                               const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                               const stationary_equation_parameters_2d<T>& auxiliary_data,
                               const bool is_neumann = false) {
    std::vector<T> right_part(mesh.container().nodes_count() + is_neumann, T{0});
    boundary_condition_second_kind_2d(right_part, mesh, boundaries_conditions);
    if (auxiliary_data.right_part)
        integrate_right_part(right_part, mesh, auxiliary_data.right_part);                                 
    if (is_neumann) {
        if (std::abs(std::reduce(right_part.begin(), right_part.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
            throw std::domain_error{"It's unsolvable Neumann problem."};
        right_part[right_part.size() - 1] = auxiliary_data.energy;
    }
    return right_part;
}

template<std::floating_point T>
std::unique_ptr<slae::preconditioner_base<T>> init_preconditioner(problem_settings settings,
                                                                  const mesh::mesh_2d<T>& mesh,
                                                                  const evaluated_conductivity_2d<T>& parameters,
                                                                  const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                  const std::vector<T>& temperature = {}) {
    const auto theroires_setter = std::views::all(mesh.container().groups_2d()) |
                                  std::views::transform([](const std::string& group) { return std::pair{group, theory_t::LOCAL}; });
    settings.theories = std::unordered_map<std::string, theory_t>(theroires_setter.begin(), theroires_setter.end());
    conductivity_matrix_2d<T> local_conductivity{mesh};
    local_conductivity.processing_nodes = std::ranges::iota_view{0zu, mesh.container().nodes_count()};
    local_conductivity.compute(parameters, settings);
    convection_condition_2d(local_conductivity.matrix(), mesh, boundaries_conditions, settings.is_inner_nodes); // TODO: fix computational range for MPI
    remove_first_kind_elements(local_conductivity.matrix(), settings.is_inner_nodes);
    if (!temperature.empty())
        radiation_condition_2d(local_conductivity.matrix(), mesh, boundaries_conditions, temperature, settings.is_inner_nodes);
    return slae::init_preconditioner(std::move(local_conductivity.matrix()), settings.is_symmetric());
}

template<std::floating_point T>
std::vector<T> stationary_heat_equation_solver_2d_linear(const problem_settings& settings,
                                                         const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                         const evaluated_conductivity_2d<T>& parameters,
                                                         const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                         std::vector<T> right_part) {
    conductivity_matrix_2d<T> conductivity{*mesh};
    conductivity.compute(parameters, settings);
    convection_condition_2d(conductivity.matrix(), *mesh, boundaries_conditions, settings.is_inner_nodes);
    if (!settings.is_neumann)
        boundary_condition_first_kind_2d(conductivity.matrix(), right_part, settings, mesh->container(), boundaries_conditions);
    auto solver = slae::init_iterative_solver(conductivity.matrix(), settings.is_symmetric());
    if (settings.is_nonlocal())
        solver->preconditioner(init_preconditioner(settings, *mesh, parameters, boundaries_conditions));
    return solver->solve(right_part);
}

template<std::floating_point T>
std::vector<T> stationary_heat_equation_solver_2d_nonlinear(const problem_settings& settings,
                                                            const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                            const parameters_2d<T>& parameters,
                                                            const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                            const stationary_equation_parameters_2d<T>& auxiliary_data,
                                                            const std::vector<T>& initial_right_part) {
    std::vector<T> temperature(mesh->container().nodes_count() + settings.is_neumann, T{0});
    if (auxiliary_data.initial_distribution)
        for(const size_t node : mesh->container().nodes())
            temperature[node] = auxiliary_data.initial_distribution(mesh->container().node_coord(node));
    first_kind_fill_2d(temperature, mesh->container(), boundaries_conditions, true);
    std::vector<T> delta_temperature = temperature;

    conductivity_matrix_2d<T> conductivity{*mesh};
    T difference = std::numeric_limits<T>::max();
    uintmax_t iteration = 0;
    do {
        using namespace metamath::operators;
        // TODO: Need to implement updating only nonlinear part of matrix and parameters logic
        const auto conductivity_parameters = evaluate_conductivity(*mesh, parameters, mesh::utils::nodes_to_qnodes<T>(*mesh, temperature));
        conductivity.compute(conductivity_parameters, settings);

        std::vector<T> right_part = initial_right_part;
        radiation_condition_2d(right_part, *mesh, boundaries_conditions, temperature, settings.is_inner_nodes);
        if (!settings.is_neumann)
            boundary_condition_first_kind_2d(conductivity.matrix(), right_part, settings, mesh->container(), boundaries_conditions);
        right_part -= settings.is_symmetric() ? 
                      conductivity.matrix().template self_adjoint<metamath::linear::matrix_part::Upper>() * temperature :
                      conductivity.matrix() * temperature;
        first_kind_fill_2d(right_part, mesh->container(), boundaries_conditions, false);

        convection_condition_2d(conductivity.matrix(), *mesh, boundaries_conditions, settings.is_inner_nodes);
        radiation_condition_2d(conductivity.matrix(), *mesh, boundaries_conditions, temperature, settings.is_inner_nodes);

        // TODO: Implement smart invalidation of solver.
        auto solver = slae::init_iterative_solver(conductivity.matrix(), settings.is_symmetric());
        if (settings.is_nonlocal())
            solver->preconditioner(init_preconditioner(settings, *mesh, conductivity_parameters, boundaries_conditions));
        delta_temperature = solver->solve(right_part, temperature);

        temperature += delta_temperature;
        difference = metamath::linear::norm(delta_temperature) / (metamath::linear::norm(temperature) ?: T{1});
        ++iteration;
    } while (iteration < auxiliary_data.max_iterations && difference > auxiliary_data.tolerance);
    return temperature;
}

template<std::floating_point T>
heat_equation_solution_2d<T> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T>>& mesh,
                                                                const parameters_2d<T>& parameters,
                                                                const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                const stationary_equation_parameters_2d<T>& auxiliary_data) {
    static constexpr bool Is_Stationary = true;
    const auto settings = init_problem_settings(mesh->container(), parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);
    std::vector<T> right_part = init_right_part(*mesh, boundaries_conditions, auxiliary_data, settings.is_neumann);
    evaluated_conductivity_2d<T> conductivity_parameters;
    std::vector<T> temperature;
    if(settings.is_nonlinear()) {
        temperature = stationary_heat_equation_solver_2d_nonlinear(settings, mesh, parameters, boundaries_conditions, auxiliary_data, right_part);
        conductivity_parameters = evaluate_conductivity(*mesh, parameters, mesh::utils::nodes_to_qnodes<T>(*mesh, temperature));
    } else {
        conductivity_parameters = evaluate_conductivity(*mesh, parameters, {});
        temperature = stationary_heat_equation_solver_2d_linear(settings, mesh, conductivity_parameters, boundaries_conditions, std::move(right_part));
    }
    auto solution = heat_equation_solution_2d<T>{mesh, conductivity_parameters, temperature};
    solution.calc_flux();
    return solution;
}

}
#pragma once

#include "save_results.hpp"

#include <config/thermal_auxiliary_data.hpp>
#include <config/time_data.hpp>
#include <logger/logger.hpp>
#include <mesh/mesh_2d/cuthill_mckee.hpp>
#include <mesh/mesh_2d/find_neighbours.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>
#include <solvers/solver_2d/thermal/nonstationary_heat_equation_solver_2d.hpp>

namespace nonlocal {

template<std::floating_point T, std::integral I>
void save_solution(const solver_2d::thermal::heat_equation_solution_2d<T, I>& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (parallel::MPI_rank() != 0) // Only the master process saves data
        return;
    if (step.has_value())
        logger::info() << "step = " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + save.get_name("csv", "solution"), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh().container(), 
        {{"temperature", solution.temperature()}, {"flux_x", solution.flux()[X]}, {"flux_y", solution.flux()[Y]}},
        save.precision()
    );
}

template<std::floating_point T, std::signed_integral I>
solver_2d::thermal::heat_equation_solution_2d<T> solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const solver_2d::thermal::parameters_2d<T>& parameters,
    const solver_2d::thermal::thermal_boundaries_conditions_2d<T>& boundaries_conditions,
    const config::thermal_auxiliary_data_2d<T>& auxiliary) {
    const solver_2d::thermal::stationary_equation_parameters_2d<T> auxiliary_data {
        .right_part = [value = std::get<spatial_dependency<T, 2>>(auxiliary.right_part)](const std::array<T, 2>& x) constexpr noexcept { return value(x); },
        .initial_distribution = [value = auxiliary.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return value; },
        .tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-10,
        .max_iterations = 40,
        .energy = auxiliary.energy
    };
    auto solution = solver_2d::thermal::stationary_heat_equation_solver_2d<I>( 
        mesh, parameters, boundaries_conditions, auxiliary_data
    );
    solution.calc_flux();
    return solution;
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const solver_2d::thermal::parameters_2d<T>& parameters,
    const solver_2d::thermal::thermal_boundaries_conditions_2d<T>& boundaries_conditions,
    const config::thermal_auxiliary_data_2d<T>& auxiliary,
    const config::time_data<T>& time,
    const config::save_data& save) {
    solver_2d::thermal::nonstationary_heat_equation_solver_2d<T, uint32_t, I> solver{mesh, time.time_step};
    solver.compute(parameters, boundaries_conditions,
        [init_dist = auxiliary.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return init_dist; });
    {
        solver_2d::thermal::heat_equation_solution_2d<T> solution{mesh, parameters, solver.temperature()};
        solution.calc_flux();
        save_solution(solution, save, 0u);
    }
    for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
        solver.calc_step(boundaries_conditions, 
        [value = std::get<spatial_dependency<T, 2>>(auxiliary.right_part)](const std::array<T, 2>& x) constexpr noexcept { return value(x); });
        if (step % time.save_frequency == 0) {
            solver_2d::thermal::heat_equation_solution_2d<T> solution{mesh, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, save, step);
        }
    }
}

}
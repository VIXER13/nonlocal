#pragma once

#include "save_results.hpp"

#include "config/thermal_auxiliary_data.hpp"
#include "config/time_data.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_2d.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"
#include "cuthill_mckee.hpp"
#include "find_neighbours.hpp"

namespace nonlocal::thermal {

template<std::floating_point T, std::integral I>
void save_solution(const heat_equation_solution_2d<T, I>& solution, 
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
heat_equation_solution_2d<T> solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const parameters_2d<T>& parameters,
    const thermal::thermal_boundaries_conditions_2d<T>& boundaries_conditions,
    const config::thermal_auxiliary_data<T>& auxiliary) {
    auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
        mesh, parameters, boundaries_conditions, 
        [value = auxiliary.right_part](const std::array<T, 2>& x) constexpr noexcept { return value; },
        auxiliary.energy
    );
    solution.calc_flux();
    return solution;
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const parameters_2d<T>& parameters,
    const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
    const config::thermal_auxiliary_data<T>& auxiliary,
    const config::time_data<T>& time,
    const config::save_data& save) {
    nonstationary_heat_equation_solver_2d<T, uint32_t, I> solver{mesh, time.time_step};
    solver.compute(parameters, boundaries_conditions,
        [init_dist = auxiliary.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return init_dist; });
    {
        nonlocal::thermal::heat_equation_solution_2d<T> solution{mesh, parameters, solver.temperature()};
        solution.calc_flux();
        save_solution(solution, save, 0u);
    }
    for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
        solver.calc_step(boundaries_conditions,
            [right_part = auxiliary.right_part](const std::array<T, 2>& x) constexpr noexcept { return right_part; });
        if (step % time.save_frequency == 0) {
            nonlocal::thermal::heat_equation_solution_2d<T> solution{mesh, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, save, step);
        }
    }
}

}
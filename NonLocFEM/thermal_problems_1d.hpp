#pragma once

#include "save_results.hpp"

#include <config/thermal_auxiliary_data.hpp>
#include <config/time_data.hpp>
#include <config/read_mesh.hpp>
#include <config/read_thermal_boundary_conditions.hpp>
#include <config/read_thermal_parameters.hpp>
#include <solvers/solver_1d/thermal/stationary_heat_equation_solver_1d.hpp>
#include <solvers/solver_1d/thermal/nonstationary_heat_equation_solver_1d.hpp>

namespace nonlocal {

template<std::floating_point T>
void save_solution(const solver_1d::thermal::heat_equation_solution_1d<T>& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (step.has_value())
        logger::info() << "save step " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + save.get_name("csv", "solution"), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh(), {{"temperature", solution.temperature()}, {"flux", solution.flux()}}, save.precision());
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_1d_problem(const nlohmann::json& config, const config::save_data& save, const bool time_dependency) {
    const auto mesh = config::read_mesh_1d<T>(config, {});
    const auto parameters = config::read_thermal_parameters_1d<T>(config["materials"], "materials");
    const auto auxiliary = config::thermal_auxiliary_data_1d<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    auto boundaries_conditions = config::read_thermal_boundaries_conditions_1d<T>(config["boundaries"], "boundaries");

    if (!time_dependency) {
        auto solution = solver_1d::thermal::stationary_heat_equation_solver_1d<T, I>(
            mesh, parameters, boundaries_conditions,
            solver_1d::thermal::stationary_equation_parameters_1d<T>{
                .right_part = [value = auxiliary.right_part](const T x) constexpr noexcept { return value; },
                .initial_distribution = [value = auxiliary.initial_distribution](const T x) constexpr noexcept { return value; },
                .energy = auxiliary.energy
            }
        );
        solution.calc_flux();
        save_solution(solution, save);
    } else {
        //const auto flux = [](const T t) { return 1048576. / 315. * std::exp(-8 * t) * metamath::functions::power<8>(t); };
        //boundaries_conditions.front() = std::make_unique<solver_1d::thermal::flux_1d<T>>(flux(0));

        config::check_required_fields(config, {"time"});
        const config::time_data<T> time{config["time"], "time"};
        solver_1d::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, parameters, time.time_step};
        solver.initialize_temperature(std::vector(mesh->nodes_count(), auxiliary.initial_distribution));
        solver.compute(boundaries_conditions);
        solver_1d::thermal::heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature(), time.time_step};
        solution.calc_flux();
        save_solution(solution, save, 0u);
        for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
            //boundaries_conditions.front() = std::make_unique<solver_1d::thermal::flux_1d<T>>(flux(step * time.time_step));
            solver.calc_step(boundaries_conditions,
                [right_part = auxiliary.right_part](const T x) constexpr noexcept { return right_part; });
            solution.temperature(solver.temperature());
            solution.time(solver.time());
            solution.calc_flux();
            if (step % time.save_frequency == 0) {
                save_solution(solution, save, step);
            }
        }
    }
}

}
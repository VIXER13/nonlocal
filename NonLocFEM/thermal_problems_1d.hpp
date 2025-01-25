#pragma once

#include "config/read_mesh_1d.hpp"
#include "config/read_thermal_boundary_conditions.hpp"

#include "problems_utils.hpp"

#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
//#include "thermal/nonstationary_relax_time_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

namespace nonlocal::thermal {

template<std::floating_point T>
parameters_1d<T> make_thermal_parameters(
    const typename config::thermal_materials_1d<T>::materials_t& materials) {
    parameters_1d<T> parameters(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        parameters[i] = {
            .model = {
                .influence = influence::polynomial_1d<T, 1, 1>{materials[i].model.nonlocal_radius},
                .local_weight = materials[i].model.local_weight
            },
            .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(
                materials[i].physical.conductivity,
                materials[i].physical.capacity,
                materials[i].physical.density,
                materials[i].physical.relaxation_time 
            )
        };
    return parameters;
}

template<std::floating_point T>
void save_solution(const thermal::heat_equation_solution_1d<T>& solution, 
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
    const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
    const auto mesh = nonlocal::config::read_mesh_1d<T>(config, {});
    const auto parameters = make_thermal_parameters<T>(materials.materials);
    const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    const auto boundaries_conditions = config::read_thermal_boundaries_conditions_1d<T>(config["boundaries"], "boundaries");

    if (!time_dependency) {
        auto solution = stationary_heat_equation_solver_1d<T, I>(
            mesh, parameters, boundaries_conditions,
            stationary_equation_parameters_1d<T>{
                .right_part = [value = auxiliary.right_part](const T x) constexpr noexcept { return value; },
                .initial_distribution = [value = auxiliary.initial_distribution](const T x) constexpr noexcept { return value; },
                .energy = auxiliary.energy
            }
        );
        solution.calc_flux();
        save_solution(solution, save);
    } else {
        config::check_required_fields(config, {"time"});
        const config::time_data<T> time{config["time"], "time"};
        //nonstationary_relax_time_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};
        nonstationary_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = auxiliary.initial_distribution](const T x) constexpr noexcept { return init_dist; });
        {   // Step 0
            heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, save, 0u);
        }
        //std::vector<T> relaxation_integral(solver.temperature().size());
        for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = auxiliary.right_part](const T x) constexpr noexcept { return right_part; });
            heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
            // if (solver._relaxation_time) {
            //     const T time = step * solver.time_step();
            //     using namespace metamath::functions;
            //     relaxation_integral *= std::exp(-solver.time_step() / solver._relaxation_time);
            //     relaxation_integral += (solver.time_step() / solver._relaxation_time) * solution.calc_flux();
            //     solution.calc_relaxation_flux(relaxation_integral, time, solver._relaxation_time);
            // }
            if (step % time.save_frequency == 0) {
                solution.calc_flux();
                save_solution(solution, save, step);
            }
        }
    }
}

}
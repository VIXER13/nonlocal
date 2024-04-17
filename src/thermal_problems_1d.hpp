#ifndef NONLOCFEM_THERMAL_PROBLEMS_1D_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_1D_HPP

#include "problems_utils.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/heat_equation_solver_with_relaxation_1d.hpp"
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
std::unique_ptr<thermal_boundary_condition_1d<T>> make_thermal_boundary_condition(
    const config::thermal_boundary_condition_data<T, 1>& condition) {
    switch (condition.kind) {
    case config::thermal_boundary_condition_t::TEMPERATURE:
        return std::make_unique<temperature_1d<T>>(condition.temperature);

    case config::thermal_boundary_condition_t::FLUX:
        return std::make_unique<flux_1d<T>>(condition.flux);

    case config::thermal_boundary_condition_t::CONVECTION:
        return std::make_unique<convection_1d<T>>(condition.heat_transfer, condition.temperature);

    case config::thermal_boundary_condition_t::RADIATION:
        return std::make_unique<radiation_1d<T>>(condition.emissivity, T{0});

    case config::thermal_boundary_condition_t::COMBINED:
        return std::make_unique<combined_flux_1d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity, T{0});

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

template<std::floating_point T>
thermal_boundaries_conditions_1d<T> make_thermal_boundaries_conditions_1d(
    const config::thermal_boundaries_conditions_1d<T>& conditions) {
    return {
        make_thermal_boundary_condition(conditions.conditions.at("left")),
        make_thermal_boundary_condition(conditions.conditions.at("right"))
    };
}

template<std::floating_point T>
void save_solution(const thermal::heat_equation_solution_1d<T>& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (step.has_value())
        logger::get().log() << "save step " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + save.get_name("csv", "solution"), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh(), {{"temperature", solution.temperature()}, {"flux", solution.flux()}}, save.precision());
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_1d_problem(const nlohmann::json& config, const config::save_data& save, const bool time_dependency) {
    const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
    const auto mesh = make_mesh_1d(get_segments_data(materials), 
                                   config::mesh_data<1u>{config.value("mesh", nlohmann::json::object()), "mesh"});
    const auto parameters = make_thermal_parameters<T>(materials.materials);
    const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    auto boundaries_conditions = make_thermal_boundaries_conditions_1d(
        config::thermal_boundaries_conditions_1d<T>{config["boundaries"], "boundaries"}
    );

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
        nonstationary_heat_equation_solver_1d<T, I> solver{mesh, parameters, time.time_step};

        //static constexpr auto exact_solution = [](const T t, const T x) { return x * t; }; // example 1
        //static constexpr auto exact_solution = [](const T t, const T x) { return x * x * t; }; // example 2
        //static constexpr auto exact_solution = [](const T t, const T x) { return x * x * t * t; }; // example 3
        //static constexpr auto exact_solution = [](const T t, const T x) { return std::exp(x) * t; }; // example 4
        //static constexpr auto exact_solution = [](const T t, const T x) { return std::exp(-x) * (t + 1); }; // example 5
        static constexpr auto exact_solution = [](const T t, const T x) { return std::exp(-t) * (x + 1) * (x + 1); }; // example 6
        
        solver.initialize_temperature([](const T x) { return exact_solution(0, x); }); // example N

        solver.compute(boundaries_conditions);
        {   // Step 0
            heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, save, 0u);
        }
        for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
            const T t = time.time_step * step; // example N
            boundaries_conditions.front() = std::make_unique<temperature_1d<T>>(exact_solution(t, 0)); // example N
            boundaries_conditions.back() = std::make_unique<temperature_1d<T>>(exact_solution(t, 1)); // example N
            //const auto right_part = [](const T) { return 0; };
            //const auto right_part = [t](const T x) { return x; }; // example 1
            //const auto right_part = [t](const T x) { return 2 - 2 * std::exp(-t) - 2 * t + x * x; }; // example 2
            //const auto right_part = [t](const T x) { return -4 + 4 * std::exp(-t) + 4 * t - 2 * t * t + 2 * t * x * x; }; // example 3
            //const auto right_part = [t](const T x) { return 2 * std::exp(x) - std::exp(-t + x) - std::exp(x) * t; }; // example 4
            //const auto right_part = [t](const T x) { return std::exp(-x) * (1 - t); }; // example 5
            const auto right_part = [t](const T x) { return -std::exp(-t) * (1 + 2 * (t + x) + x * x); }; // example 6

            solver.calc_step(boundaries_conditions, right_part);
            if (step % time.save_frequency == 0) {
                heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
                solution.calc_flux();
                save_solution(solution, save, step);
            }
        }
    }
}

}

#endif
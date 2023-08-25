#ifndef NONLOCFEM_THERMAL_PROBLEMS_UTILS_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_UTILS_HPP

#include "problems_utils.hpp"

#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

namespace nonlocal::thermal {

template<std::floating_point T>
parameters_1d<T> make_thermal_parameters(
    const typename config::thermal_materials_1d<T>::materials_t& materials) {
    parameters_1d<T> parameters(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        parameters[i] = {
            .model = {
                .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{materials[i].model.nonlocal_radius},
                .local_weight = materials[i].model.local_weight
            },
            .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(
                materials[i].physical.conductivity,
                materials[i].physical.capacity,
                materials[i].physical.density
            )
        };
    return parameters;
}

template<std::floating_point T>
std::unique_ptr<thermal_boundary_condition_1d<T>> make_thermal_boundary_condition(
    const config::thermal_boundary_condition_data<T, 1>& condition) {
    switch (condition.kind) {
    case nonlocal::config::thermal_boundary_condition_t::TEMPERATURE:
        return std::make_unique<temperature_1d<T>>(condition.temperature);

    case nonlocal::config::thermal_boundary_condition_t::FLUX:
        return std::make_unique<flux_1d<T>>(condition.flux);

    case nonlocal::config::thermal_boundary_condition_t::CONVECTION:
        return std::make_unique<convection_1d<T>>(condition.heat_transfer, condition.temperature);

    case nonlocal::config::thermal_boundary_condition_t::RADIATION:
        return std::make_unique<radiation_1d<T>>(condition.emissivity, T{0});

    case nonlocal::config::thermal_boundary_condition_t::COMBINED:
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
void save_solution(thermal::heat_equation_solution_1d<T>&& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (step);
        std::cout << "INFO: step = " << *step << std::endl;
    const auto save_vector = [&solution, &save, step](const std::vector<T>& x, const std::string& name) {
        const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + name, "csv") : save.path(name, "csv");
        nonlocal::mesh::utils::save_as_csv(solution.mesh(), x, path, save.precision());
    };
    if (save.contains("temperature"))
        save_vector(solution.temperature(), "temperature");
    if (save.contains("flux"))
        save_vector(solution.calc_flux(), "flux");
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_1d_problem(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem) {
    const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
        const auto mesh = make_mesh_1d(get_segments_data(materials), 
                                       config::mesh_data<1u>{config.value("mesh", nlohmann::json::object()), "mesh"});
        const auto parameters = thermal::make_thermal_parameters<T>(materials.materials);
        const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
        const auto boundaries_conditions = thermal::make_thermal_boundaries_conditions_1d(
            config::thermal_boundaries_conditions_1d<T>{config["boundaries"], "boundaries"}
        );

        if (problem == config::problem_t::THERMAL_STATIONARY) {
            auto solution = thermal::stationary_heat_equation_solver_1d<T, I>(
                mesh, parameters, boundaries_conditions,
                thermal::stationary_equation_parameters_1d<T>{
                    .right_part = [value = auxiliary.right_part](const T x) constexpr noexcept { return value; },
                    .initial_distribution = [value = auxiliary.initial_distribution](const T x) constexpr noexcept { return value; },
                    .energy = auxiliary.energy
                }
            );
            save_solution(std::move(solution), save);
        } else if (problem == config::problem_t::THERMAL_NONSTATIONARY) {
            config::check_required_fields(config, {"time"});
            const config::time_data<T> time{config["time"], "time"};
            thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};
            solver.compute(parameters, boundaries_conditions,
                [init_dist = auxiliary.initial_distribution](const T x) constexpr noexcept { return init_dist; });
            save_solution(thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, save, 0u);
            for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
                solver.calc_step(boundaries_conditions,
                    [right_part = auxiliary.right_part](const T x) constexpr noexcept { return right_part; });
                if (step % time.save_frequency == 0)
                    save_solution(thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, save, step);
            }
        }
}

}

#endif
#ifndef NONLOCFEM_THERMAL_PROBLEMS_2D_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_2D_HPP

#include "problems_utils.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_2d.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"
#include "influence_functions_2d.hpp"

namespace nonlocal::thermal {

template<std::floating_point T>
parameters_2d<T> make_parameters(const config::thermal_materials_2d<T>& materials) {
    parameters_2d<T> parameters;
    for(const auto& [name, material] : materials.materials) {
        auto& parameter = parameters[name] = {
            .model = {
                .influence = influence::polynomial_2d<T, 2, 1>{material.model.nonlocal_radius},
                .local_weight = material.model.local_weight
            },
            .physical = std::make_shared<parameter_2d<T, coefficients_t::CONSTANTS>>(
                metamath::types::make_square_matrix<T, 2u>(material.physical.conductivity),
                material.physical.capacity,
                material.physical.density
            )
        };
        parameter.physical->material = material.physical.material;
    }
    return parameters;
}

template<std::floating_point T>
thermal_boundary_condition_2d<T> make_boundary_condition(const config::thermal_boundary_condition_data<T, 2>& condition) {
    switch (condition.kind) {
    case config::thermal_boundary_condition_t::TEMPERATURE:
        return std::make_unique<temperature_2d<T>>(condition.temperature);

    case config::thermal_boundary_condition_t::FLUX:
        return std::make_unique<flux_2d<T>>(condition.flux);

    case config::thermal_boundary_condition_t::CONVECTION:
        return std::make_unique<convection_2d<T>>(condition.heat_transfer, condition.temperature);

    case config::thermal_boundary_condition_t::RADIATION:
        return std::make_unique<radiation_2d<T>>(condition.emissivity);

    case config::thermal_boundary_condition_t::COMBINED:
        return std::make_unique<combined_flux_2d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity);

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

template<std::floating_point T>
thermal_boundaries_conditions_2d<T> make_boundaries_conditions(const config::thermal_boundaries_conditions_2d<T>& conditions) {
    thermal_boundaries_conditions_2d<T> result;
    for(const auto& [name, condition] : conditions.conditions)
        result[name] = make_boundary_condition(condition);
    return result;
}

template<std::floating_point T, std::signed_integral I>
void save_solution(heat_equation_solution_2d<T, I>&& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (step);
        logger::get().log(logger::log_level::INFO) << "step = " << *step << std::endl;
    const auto save_vector = [&solution, &save, step](const std::vector<T>& x, const std::string& name) {
        const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + name, "csv") : save.path(name, "csv");
        mesh::utils::save_as_csv(path, solution.mesh().container(), x, save.precision());
    };
    if (save.contains("temperature"))
        save_vector(solution.temperature(), "temperature");
    if (const bool contains_x = save.contains("flux_x"), contains_y = save.contains("flux_y"); contains_x || contains_y) {
        solution.calc_flux();
        if (contains_x)
            save_vector(solution.flux()[X], "flux_x");
        if (contains_y)
            save_vector(solution.flux()[Y], "flux_y");
    }
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const nlohmann::json& config, 
    const config::save_data& save, const bool time_dependency) {
    const config::thermal_materials_2d<T> materials{config["materials"], "materials"};
    mesh->find_neighbours(get_search_radii(materials));
    const auto parameters = make_parameters(materials);
    const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    const auto boundaries_conditions = make_boundaries_conditions(
        config::thermal_boundaries_conditions_2d<T>{config["boundaries"], "boundaries"});
    if (!time_dependency) {
        auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
            mesh, parameters, boundaries_conditions, 
            [value = auxiliary.right_part](const std::array<T, 2>& x) constexpr noexcept { return value; },
            auxiliary.energy
        );
        save_solution(std::move(solution), save);
    } else {
        const config::time_data<T> time{config["time"], "time"};
        nonstationary_heat_equation_solver_2d<T, I, int64_t> solver{mesh, time.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = auxiliary.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return init_dist; });
        save_solution(nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, parameters, solver.temperature()}, save, 0u);
        for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = auxiliary.right_part](const std::array<T, 2>& x) constexpr noexcept { return right_part; });
            if (step % time.save_frequency == 0)
                save_solution(nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, parameters, solver.temperature()}, save, step);
        }
    }
}

}

#endif
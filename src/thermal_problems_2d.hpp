#ifndef NONLOCFEM_THERMAL_PROBLEMS_2D_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_2D_HPP

#include "problems_utils.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_2d.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"


namespace nonlocal::thermal {

template<std::floating_point T>
parameters_2d<T> make_parameters(const config::thermal_materials_2d<T>& materials) {
    parameters_2d<T> parameters;
    for(const auto& [name, material] : materials.materials) {
        auto& parameter = parameters[name] = {
            .model = {
                .influence = get_influence(material.model.influence, material.model.nonlocal_radius),
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
void save_solution(const heat_equation_solution_2d<T, I>& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (parallel_utils::MPI_rank() != 0) // Only the master process saves data
        return;
    if (step);
        logger::get().log(logger::log_level::INFO) << "step = " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + save.get_name("csv", "solution"), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh().container(), 
        {{"temperature", solution.temperature()}, {"flux_x", solution.flux()[X]}, {"flux_y", solution.flux()[Y]}},
        save.precision()
    );
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const nlohmann::json& config, 
    const config::save_data& save, const bool time_dependency) {
    const config::thermal_materials_2d<T> materials{config["materials"], "materials"};
    mesh->find_neighbours(get_search_radii(materials));
    mesh->balancing(mesh::balancing_t::MEMORY, true);
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
        solution.calc_flux();
        save_solution(solution, save);
    } else {
        config::check_required_fields(config, {"time"});
        const config::time_data<T> time{config["time"], "time"};
        nonstationary_heat_equation_solver_2d<T, I, int64_t> solver{mesh, time.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = auxiliary.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return init_dist; });
        {
            nonlocal::thermal::heat_equation_solution_2d<T, I> solution{mesh, parameters, solver.temperature()};
            solution.calc_flux();
            save_solution(solution, save, 0u);
        }
        for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = auxiliary.right_part](const std::array<T, 2>& x) constexpr noexcept { return right_part; });
            if (step % time.save_frequency == 0) {
                nonlocal::thermal::heat_equation_solution_2d<T, I> solution{mesh, parameters, solver.temperature()};
                solution.calc_flux();
                save_solution(solution, save, step);
            }
        }
    }
}

}

#endif
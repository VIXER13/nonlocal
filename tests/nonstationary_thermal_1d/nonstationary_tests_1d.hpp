#ifndef UNIT_TESTS_BOUNDARY_CONDITIONS_HPP
#define UNIT_TESTS_BOUNDARY_CONDITIONS_HPP

#include <boost/ut.hpp>
#include <tuple>
#include <filesystem>

#include "problems_utils.hpp"

#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

namespace nonstat_1d_tests {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

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
    if (step.has_value());
        logger::get().log() << "save step " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + save.get_name("csv", "solution"), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh(), {{"temperature", solution.temperature()}, {"flux", solution.flux()}}, save.precision());
}

template<std::floating_point T, std::signed_integral I>
void check_solution(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const heat_equation_solution_1d<T>& solution, T time_layer, I step, 
                    std::function<T(T, T)> ref, T eps = epsilon) {
    auto& sol = solution.temperature();
    for (std::size_t k = 0; k < sol.size(); ++k) {
        expect(lt(std::fabs(sol[k] - ref(time_layer, mesh->node_coord(k))), eps * step));
        std::cout << "step: " <<  step << std::endl;
    }
}

template<std::floating_point T, std::signed_integral I>
void save_and_calc_flux(const config::time_data<T>& time, I step, heat_equation_solution_1d<T>& solution) {
    if (step % time.save_frequency == 0) {
        solution.calc_flux();
        //save_solution(solution, nonlocal::config::save_data{}, step);
    }
}

template<std::floating_point T, std::signed_integral I>
void solve_nonstationary_thermal_1d_problem_json(const nlohmann::json& config, std::function<T(T, T)> ref_sol)  {
    const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
    const auto mesh = make_mesh_1d(get_segments_data(materials), 
                                   config::mesh_data<1u>{config.value("mesh", nlohmann::json::object()), "mesh"});
    const auto parameters = make_thermal_parameters<T>(materials.materials);
    const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    const auto boundaries_conditions = make_thermal_boundaries_conditions_1d(
        config::thermal_boundaries_conditions_1d<T>{config["boundaries"], "boundaries"}
    );

    config::check_required_fields(config, {"time"});
    const config::time_data<T> time{config["time"], "time"};

    std::function<T(T)> init_dist = [init_dist = auxiliary.initial_distribution](const T x) constexpr noexcept { return init_dist; };
    std::function<T(T, T)> right_part = [right_part = auxiliary.right_part](const T t, const T x) constexpr noexcept { return right_part; };

    nonstationary_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};
    solver.compute(parameters, boundaries_conditions, init_dist);  

    {   // Step initial
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, I(0), solution);
    }

    for(const I step : std::ranges::iota_view{1u, time.steps_count}) {
        std::function<T(T)> rp = [&right_part, &step, &time](const T x) { return right_part((step - 1) * time.time_step, x); };
        solver.calc_step(boundaries_conditions, rp);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, step, solution);
    }

    { // Step final
        std::function<T(T)> rp = [&right_part, &time](const T x) { return right_part(time.steps_count * time.time_step, x); };
        solver.calc_step(boundaries_conditions, rp);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, time.steps_count + 1, solution);
        check_solution<T, I>(mesh, solution, time.steps_count * time.time_step, time.steps_count + 1, ref_sol);
    }
}
    
}

#endif

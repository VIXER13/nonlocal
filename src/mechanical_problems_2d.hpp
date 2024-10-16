#pragma once

#include "logger.hpp"
#include "nonlocal_config.hpp"
#include "equilibrium_equation_2d.hpp"
#include "influence_functions_2d.hpp"

namespace nonlocal::mechanical {

template<std::floating_point T>
mechanical_parameters_2d<T> make_parameters(const config::mechanical_materials_2d<T>& materials) {
    mechanical_parameters_2d<T> parameters;
    for(const auto& [name, material] : materials.materials)
        parameters.materials[name] = {
            .model = {
                .influence = get_influence(material.model.influence, material.model.nonlocal_radius),
                .local_weight = material.model.local_weight
            },
            .physical = {
                material.physical.youngs_modulus,
                material.physical.poissons_ratio,
                material.physical.thermal_expansion
            }
        };
    return parameters;
}

template<std::floating_point T>
std::unique_ptr<mechanical_boundary_condition_2d<T>> make_boundary_condition(
    const std::optional<config::mechanical_boundary_condition_data<T>>& condition) {
    if (!condition)
        return nullptr;

    switch (condition->kind) {
    case config::mechanical_boundary_condition_t::DISPLACEMENT:
        return std::make_unique<displacement_2d<T>>(condition->value);
    
    case config::mechanical_boundary_condition_t::PRESSURE:
        return std::make_unique<pressure_2d<T>>(condition->value);
    
    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition->kind))};
    }
}

template<std::floating_point T>
mechanical_boundary_conditions_2d<T> make_boundary_conditions(const config::mechanical_boundary_conditions_data<T, 2u>& conditions) {
    mechanical_boundary_conditions_2d<T> result;
    result[0] = make_boundary_condition(conditions.conditions[0]);
    result[1] = make_boundary_condition(conditions.conditions[1]);
    return result;
}

template<std::floating_point T>
mechanical_boundaries_conditions_2d<T> make_boundaries_conditions(const config::mechanical_boundaries_conditions_2d<T>& conditions) {
    mechanical_boundaries_conditions_2d<T> result;
    for(const auto& [name, conditions] : conditions.conditions)
        result[name] = make_boundary_conditions(conditions);
    return result;
}

template<std::floating_point T, std::signed_integral I>
mechanical::mechanical_solution_2d<T> solve_mechanical_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const config::mechanical_materials_2d<T>& materials,
    const config::mechanical_boundaries_conditions_2d<T>& conditions,
    const std::vector<T>& delta_temperature = {}) {
    auto parameters = make_parameters(materials);
    parameters.delta_temperature = delta_temperature;
    const auto boundaries_conditions = make_boundaries_conditions(conditions);
    auto solution = nonlocal::mechanical::equilibrium_equation<I>(
        mesh, parameters, boundaries_conditions,
        [](const std::array<T, 2>&) constexpr noexcept { return std::array<T, 2>{}; }
    );
    solution.calc_strain_and_stress();
    return solution;
}
    
}
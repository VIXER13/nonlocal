#ifndef NONLOCFEM_THERMAL_PROBLEMS_UTILS_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_UTILS_HPP

#include "problems_utils.hpp"

#include "thermal/stationary_heat_equation_solver_1d.hpp"
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

}

#endif
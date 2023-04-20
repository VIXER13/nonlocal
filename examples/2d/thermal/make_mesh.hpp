#ifndef NONLOCAL_MAKE_MESH_2D_HPP
#define NONLOCAL_MAKE_MESH_2D_HPP

#include "stationary_heat_equation_solver_2d.hpp"
#include "influence_functions_2d.hpp"
#include "thermal_config_data.hpp"

namespace nonlocal {

template<class T, class I>
auto make_mesh(const std::filesystem::path& path,
               const typename config::stationary_thermal_data<T, 2>::materials_t& materials) {
    auto mesh = std::make_shared<mesh::mesh_2d<T, I>>(path);
    std::unordered_map<std::string, T> radii;
    for(const auto& [name, material] : materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            radii.emplace(name, std::max(material.model.search_radius[0],
                                         material.model.search_radius[1]));
    mesh->find_neighbours(radii);
    return mesh;
}

template<class T>
thermal::parameters_2d<T> make_parameters(const typename config::stationary_thermal_data<T, 2>::materials_t& materials) {
    thermal::parameters_2d<T> parameters;
    for(const auto& [name, material] : materials) {
        parameters[name] = {
            .model = {
                .influence = influence::polynomial_2d<T, 2, 1>{material.model.nonlocal_radius},
                .local_weight = material.model.local_weight
            },
            .physical = std::make_shared<thermal::parameter_2d<T, coefficients_t::CONSTANTS>>(
                metamath::types::make_square_matrix<T, 2u>(material.physical.conductivity),
                material.physical.capacity,
                material.physical.density,
                material.physical.material
            )
        };
    }
    return parameters;
}

template<std::floating_point T>
thermal::thermal_boundary_condition_2d<T> make_boundary_condition(
    const config::thermal_boundary_condition_data<T, 2>& condition) {
    switch (condition.kind) {
    case thermal::boundary_condition_t::TEMPERATURE:
        return std::make_unique<thermal::temperature_2d<T>>(condition.temperature);

    case thermal::boundary_condition_t::FLUX:
        return std::make_unique<thermal::flux_2d<T>>(condition.flux);

    case thermal::boundary_condition_t::CONVECTION:
        return std::make_unique<thermal::convection_2d<T>>(condition.heat_transfer, condition.temperature);

    case thermal::boundary_condition_t::RADIATION:
        return std::make_unique<thermal::radiation_2d<T>>(condition.emissivity);

    case thermal::boundary_condition_t::COMBINED:
        return std::make_unique<thermal::combined_flux_2d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity);

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_2d<T> make_boundaries_conditions(
    const config::thermal_boundaries_conditions_data<T, 2>& conditions) {
    thermal::thermal_boundaries_conditions_2d<T> result;
    for(const auto& [name, condition] : conditions.conditions)
        result[name] = make_boundary_condition(condition);
    return result;
}

}

#endif
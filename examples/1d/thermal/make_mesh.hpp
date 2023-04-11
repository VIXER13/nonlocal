#ifndef NONLOCAL_MAKE_MESH_1D_HPP
#define NONLOCAL_MAKE_MESH_1D_HPP

#include "make_element.hpp"

#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "mesh_1d_utils.hpp"
#include "thermal_config_data.hpp"

namespace nonlocal {

template<class T, template<class, size_t> class Physics>
std::shared_ptr<nonlocal::mesh::mesh_1d<T>> make_mesh(
    const std::vector<nonlocal::config::material_data<Physics, T, 1>>& materials,
    const size_t element_order, const size_t quadrature_order) {
    std::vector<nonlocal::mesh::segment_data<T>> segments(materials.size());
    std::vector<T> search_radii(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, segments.size()}) {
        segments[i] = nonlocal::mesh::segment_data<T>{
            .length = materials[i].length,
            .elements = materials[i].elements_count
        };
        search_radii[i] = materials[i].model.search_radius;
    }
    auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
        nonlocal::make_element<T>(nonlocal::element_1d_order_t(1), 
                                  nonlocal::quadrature_1d_order_t(1)), segments);
    mesh->find_neighbours(search_radii);
    return mesh;
}

template<class T>
nonlocal::thermal::parameters_1d<T> make_thermal_parameters(
    const std::vector<nonlocal::config::material_data<nonlocal::config::thermal_material_data, T, 1>>& materials) {
    nonlocal::thermal::parameters_1d<T> parameters(materials.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.size()})
        parameters[i] = {
            .model = {
                .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{materials[i].model.nonlocal_radius},
                .local_weight = materials[i].model.local_weight
            },
            .physical = std::make_shared<thermal::parameter_1d<T, coefficients_t::CONSTANTS>>(
                materials[i].physical.conductivity,
                materials[i].physical.capacity,
                materials[i].physical.density
            ),
        };
    return parameters;
}

template<std::floating_point T>
std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>> make_boundary_condition(
    const nonlocal::config::thermal_boundary_condition_data<T, 1u>& condition) {
    switch (condition.kind) {
    case nonlocal::thermal::boundary_condition_t::TEMPERATURE:
        return std::make_unique<nonlocal::thermal::temperature_1d<T>>(condition.temperature);

    case nonlocal::thermal::boundary_condition_t::FLUX:
        return std::make_unique<nonlocal::thermal::flux_1d<T>>(condition.flux);

    case nonlocal::thermal::boundary_condition_t::CONVECTION:
        return std::make_unique<nonlocal::thermal::convection_1d<T>>(condition.heat_transfer, condition.temperature);

    case nonlocal::thermal::boundary_condition_t::RADIATION:
        return std::make_unique<nonlocal::thermal::radiation_1d<T>>(condition.emissivity, T{0});

    case nonlocal::thermal::boundary_condition_t::COMBINED:
        return std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity, T{0});

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

}

#endif
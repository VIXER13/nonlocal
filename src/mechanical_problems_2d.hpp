#ifndef NONLOCFEM_MECHANICAL_PROBLEMS_2D_HPP
#define NONLOCFEM_MECHANICAL_PROBLEMS_2D_HPP

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
                .influence = influence::polynomial_2d<T, 2, 1>{material.model.nonlocal_radius},
                .local_weight = material.model.local_weight
            },
            .physical = {
                material.physical.youngs_modulus,
                material.physical.poissons_ratio
            }
        };
    return parameters;
}

// template<std::floating_point T>
// mechanical_boundaries_conditions_2d<T> make_boundary_condition(const config::thermal_boundary_condition_data<T, 2>& condition) {
//     switch (condition.kind) {
//     case config::thermal_boundary_condition_t::TEMPERATURE:
//         return std::make_unique<temperature_2d<T>>(condition.temperature);

//     case config::thermal_boundary_condition_t::FLUX:
//         return std::make_unique<flux_2d<T>>(condition.flux);

//     case config::thermal_boundary_condition_t::CONVECTION:
//         return std::make_unique<convection_2d<T>>(condition.heat_transfer, condition.temperature);

//     case config::thermal_boundary_condition_t::RADIATION:
//         return std::make_unique<radiation_2d<T>>(condition.emissivity);

//     case config::thermal_boundary_condition_t::COMBINED:
//         return std::make_unique<combined_flux_2d<T>>(
//             condition.flux,
//             condition.heat_transfer, condition.temperature,
//             condition.emissivity);

//     default:
//         throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
//     }
// }

template<std::floating_point T, std::signed_integral I>
void solve_mechanical_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const nlohmann::json& config, 
    const config::save_data& save, const config::problem_t problem) {
    const config::mechanical_materials_2d<T> materials{config["materials"], "materials"};
    mesh->find_neighbours(get_search_radii(materials));
    const auto parameters = make_parameters(materials);
}
    
}

#endif
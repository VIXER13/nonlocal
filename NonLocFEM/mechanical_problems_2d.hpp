#pragma once

#include "logger.hpp"
#include "nonlocal_config.hpp"
#include "equilibrium_equation_2d.hpp"
#include "influence_functions_2d.hpp"

namespace nonlocal::mechanical {

template<std::floating_point T>
mechanical_parameters_2d<T> make_parameters(const config::mechanical_materials_2d<T>& materials) {
    mechanical_parameters_2d<T> parameters;
    for(const auto& [name, material] : materials.materials) {
        parameters.materials[name] = {
            .model = {
                .influence = get_influence(material.model.influence, material.model.nonlocal_radius),
                .local_weight = material.model.local_weight
            },
            .physical = {
                material.physical.youngs_modulus,
                material.physical.poissons_ratio,
                material.physical.shear_modulus,
                material.physical.thermal_expansion
            }
        };
    }
    return parameters;
}

template<std::floating_point T, std::signed_integral I>
mechanical::mechanical_solution_2d<T> solve_mechanical_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const config::mechanical_materials_2d<T>& materials,
    const mechanical::mechanical_boundaries_conditions_2d<T>& boundaries_conditions,
    const std::vector<T>& delta_temperature = {}) {
    auto parameters = make_parameters(materials);
    parameters.delta_temperature = delta_temperature;
    auto solution = nonlocal::mechanical::equilibrium_equation<I>(
        mesh, parameters, boundaries_conditions,
        [](const std::array<T, 2>&) constexpr noexcept { return std::array<T, 2>{}; }
    );
    solution.calc_strain_and_stress();
    return solution;
}
    
}
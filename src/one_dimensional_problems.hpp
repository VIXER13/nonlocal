#ifndef NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP
#define NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP

#include "init_utils.hpp"
#include "make_element_1d.hpp"
#include "mesh_1d.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"

namespace nonlocal {

template<class T, template<class, size_t> class Physics>
std::vector<mesh::segment_data<T>> get_segments_data(const config::materials_data<Physics, T, 1>& material_data) {
    std::vector<mesh::segment_data<T>> segments(material_data.materials.size());
    for(const size_t i : std::ranges::iota_view{0u, segments.size()})
        segments[i] = mesh::segment_data<T>{
            .length = material_data.materials[i].length,
            .search_radius = material_data.materials[i].model.search_radius,
            .elements = material_data.materials[i].elements_count
        };
    return segments;
}

template<class T>
std::shared_ptr<mesh::mesh_1d<T>> make_mesh_1d(
    const std::vector<mesh::segment_data<T>>& segments,
    const config::mesh_data<1u>& mesh_data) {
    return std::make_shared<mesh::mesh_1d<T>>(
        make_element<T>(mesh_data.element_order, mesh_data.quadrature_order),
        segments
    );
}

template<std::floating_point T>
std::unique_ptr<thermal::thermal_boundary_condition_1d<T>> make_thermal_boundary_condition(
    const config::thermal_boundary_condition_data<T, 1>& condition) {
    switch (condition.kind) {
    case nonlocal::config::thermal_boundary_condition_t::TEMPERATURE:
        return std::make_unique<thermal::temperature_1d<T>>(condition.temperature);

    case nonlocal::config::thermal_boundary_condition_t::FLUX:
        return std::make_unique<thermal::flux_1d<T>>(condition.flux);

    case nonlocal::config::thermal_boundary_condition_t::CONVECTION:
        return std::make_unique<thermal::convection_1d<T>>(condition.heat_transfer, condition.temperature);

    case nonlocal::config::thermal_boundary_condition_t::RADIATION:
        return std::make_unique<thermal::radiation_1d<T>>(condition.emissivity, T{0});

    case nonlocal::config::thermal_boundary_condition_t::COMBINED:
        return std::make_unique<thermal::combined_flux_1d<T>>(
            condition.flux,
            condition.heat_transfer, condition.temperature,
            condition.emissivity, T{0});

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(condition.kind))};
    }
}

template<class T, class I>
void one_dimensional_problems(const nlohmann::json& config, const nonlocal::config::save_data& save, const config::problem_t problem) {
    static const std::set<config::problem_t> available_problems = {
        config::problem_t::THERMAL_STATIONARY,
        config::problem_t::THERMAL_NONSTATIONARY
    };
    if (!available_problems.contains(problem)) {
        throw std::domain_error{
            "In the one-dimensional case, the following problems are available: " +
            init_available_problems_list(available_problems)
        };
    }

    config::check_required_fields(config, {"materials", "boundaries"});
    if (is_thermal_problem(problem)) {
        const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
        const auto mesh = make_mesh_1d(get_segments_data(materials), 
                                       config::mesh_data<1u>{config.value("mesh", nlohmann::json::object())});
        const config::thermal_boundaries_conditions_1d<T> conditions{config["boundaries"]};
        const auto left = make_thermal_boundary_condition<T>(conditions.conditions.at("left"));

        if (problem == config::problem_t::THERMAL_STATIONARY) {
        
            std::cout << "thermal_stationary_1d" << std::endl;
        } else if (problem == config::problem_t::THERMAL_NONSTATIONARY)
            std::cout << "thermal_nonstationary_1d" << std::endl;
    }
    
}
    
}

#endif
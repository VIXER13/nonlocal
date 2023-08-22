#ifndef NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP
#define NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP

#include "init_utils.hpp"

namespace nonlocal {

/*
template<class T, template<class, size_t> class Physics>
std::shared_ptr<nonlocal::mesh::mesh_1d<T>> make_mesh_1d(
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
        nonlocal::make_element<T>(nonlocal::element_1d_order_t(element_order), 
                                  nonlocal::quadrature_1d_order_t(quadrature_order)), segments);
    mesh->find_neighbours(search_radii);
    return mesh;
}
*/

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

    const config::mesh_data<1u> mesh_info{config.value("mesh", nlohmann::json::object())};
    if (problem == config::problem_t::THERMAL_STATIONARY) {
        config::check_required_fields(config, {"materials"});
        const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
        std::cout << "thermal_stationary_1d" << std::endl;
    } else if (problem == config::problem_t::THERMAL_NONSTATIONARY)
        std::cout << "thermal_nonstationary_1d" << std::endl;
}
    
}

#endif
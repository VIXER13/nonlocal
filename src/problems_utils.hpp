#ifndef NONLOCFEM_PROBLEMS_UTILS_HPP
#define NONLOCFEM_PROBLEMS_UTILS_HPP

#include "make_element_1d.hpp"

#include "nonlocal_config.hpp"
#include "mesh_1d.hpp"

#include <set>

namespace nonlocal {

void init_save_data(const nonlocal::config::save_data& save, const nlohmann::json& config);
std::string init_available_problems_list(const std::set<config::problem_t>& available_problems);

bool is_thermal_problem(const config::problem_t problem);
bool is_mechanical_problem(const config::problem_t problem);

template<std::floating_point T, template<class, size_t> class Physics>
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

template<std::floating_point T>
std::shared_ptr<mesh::mesh_1d<T>> make_mesh_1d(
    const std::vector<mesh::segment_data<T>>& segments,
    const config::mesh_data<1u>& mesh_data) {
    return std::make_shared<mesh::mesh_1d<T>>(
        make_element<T>(mesh_data.element_order, mesh_data.quadrature_order),
        segments
    );
}

}

#endif
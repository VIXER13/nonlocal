#ifndef TESTS_PROBLEMS_UTILS_1D_HPP
#define TESTS_PROBLEMS_UTILS_1D_HPP

#include "make_element_1d.hpp"

#include "nonlocal_config.hpp"
#include "mesh_1d.hpp"


namespace nonstat_1d_tests {

using namespace nonlocal;
using namespace nonlocal::thermal;

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

template<std::floating_point T>
std::shared_ptr<mesh::mesh_1d<T>> make_mesh_1d(
    const std::vector<mesh::segment_data<T>>& segments,
    const config::order_t element_order, const config::order_t& quadrature_order) {
    return std::make_shared<mesh::mesh_1d<T>>(
        make_element<T>(element_order, quadrature_order),
        segments
    );
}

template<std::floating_point T, template<class, size_t> class Physics>
std::unordered_map<std::string, T> get_search_radii(const config::materials_data<Physics, T, 2>& materials) {
    std::unordered_map<std::string, T> result;
    for(const auto& [name, material] : materials.materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            result[name] = std::max(material.model.search_radius[X], material.model.search_radius[Y]);
    return result;
}

}

#endif
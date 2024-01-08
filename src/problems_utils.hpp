#ifndef NONLOCFEM_PROBLEMS_UTILS_2D_HPP
#define NONLOCFEM_PROBLEMS_UTILS_2D_HPP

#include "make_element_1d.hpp"

#include "influence_functions_2d.hpp"
#include "nonlocal_config.hpp"
#include "mesh_1d.hpp"

namespace nonlocal {

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

template<std::floating_point T, template<class, size_t> class Physics>
std::unordered_map<std::string, T> get_search_radii(const config::materials_data<Physics, T, 2>& materials) {
    std::unordered_map<std::string, T> result;
    for(const auto& [name, material] : materials.materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            result[name] = std::max(material.model.search_radius[X], material.model.search_radius[Y]);
    return result;
}

template<std::floating_point T>
std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> get_influence(
    const config::influence_data<T, 2>& data, const std::array<T, 2>& r) {
    if (data.family == config::influence_family_t::CUSTOM)
        throw std::domain_error{"Сustom influence functions are not currently supported."};
    if (data.family == config::influence_family_t::CONSTANT) {
        if (data.n == std::numeric_limits<size_t>::max())
            return influence::constant_2d<T, std::numeric_limits<size_t>::max()>{r};
        return influence::constant_2d_dynamic<T>{r, T(data.n)};
    }
    if (data.family == config::influence_family_t::POLYNOMIAL) {
        if (data.n == 2 && data.p == 2 && data.q == 1)
            return influence::polynomial_2d<T, 2, 1>{r};
        return influence::polynomial_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
    }
    if (data.family == config::influence_family_t::EXPONENTIAL) {
        if (data.n == 2 && data.p == 2 && data.q == 0) // q == 0 temporary condition. In future refactoring will be changed 
            return influence::normal_distribution_2d<T>{r};
        return influence::exponential_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
    }
    throw std::domain_error{"Unsupported influence family."};
}

}

#endif
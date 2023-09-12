#ifndef NONLOCFEM_THERMAL_PROBLEMS_2D_HPP
#define NONLOCFEM_THERMAL_PROBLEMS_2D_HPP

#include "nonlocal_config.hpp"
#include "thermal/stationary_heat_equation_solver_2d.hpp"
#include "influence_functions_2d.hpp"

namespace nonlocal::thermal {

template<std::floating_point T, template<class, size_t> class Physics>
std::unordered_map<std::string, T> get_search_radii(const config::materials_data<Physics, T, 2>& materials) {
    std::unordered_map<std::string, T> result;
    for(const auto& [name, material] : materials.materials)
        if (theory_type(material.model.local_weight) == theory_t::NONLOCAL)
            result[name] = std::max(material.model.search_radius[X], material.model.search_radius[Y]);
    return result;
}

template<std::floating_point T, std::signed_integral I>
auto make_mesh(const std::filesystem::path& path, const std::unordered_map<std::string, T>& search_radii) {
    auto mesh = std::make_shared<mesh::mesh_2d<T, I>>(path);
    mesh->find_neighbours(search_radii);
    return mesh;
}

template<std::floating_point T, std::signed_integral I>
void solve_thermal_2d_problem(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem) {
    const config::mesh_data<2> mesh_data{config["mesh"], "mesh"};
    const config::thermal_materials_2d<T> materials{config["materials"], "materials"};
    const auto mesh = make_mesh<T, I>(mesh_data.path, get_search_radii(materials));
}

}

#endif
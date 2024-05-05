#pragma once

#include "make_element_1d.hpp"

#include "influence_functions_2d.hpp"
#include "nonlocal_config.hpp"
#include "mesh_1d.hpp"
#include "thermal/heat_equation_solution_2d.hpp"
#include "mechanical/mechanical_solution_2d.hpp"

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
        if (data.n == 1)
            return influence::constant_2d<T, 1>{r};
        if (data.n == 2)
            return influence::constant_2d<T, 2>{r};
        if (data.n == std::numeric_limits<size_t>::max())
            return influence::constant_2d<T, std::numeric_limits<size_t>::max()>{r};
        return influence::constant_2d_dynamic<T>{r, T(data.n)};
    }
    if (data.family == config::influence_family_t::POLYNOMIAL) {
        if (data.n == 2 && data.p == 2 && data.q == 1)
            return influence::polynomial_2d<T, 2, 1>{r};
        if (data.n == std::numeric_limits<size_t>::max() && data.p == 1 && data.q == 1)
            return influence::polynomial_2d<T, 1, 1, std::numeric_limits<size_t>::max()>{r};
        return influence::polynomial_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
    }
    if (data.family == config::influence_family_t::EXPONENTIAL) {
        // q == 0 temporary condition. In future refactoring will be changed 
        if (data.n == 1 && data.p == 2 && data.q == 0)
            return influence::exponential_2d<T, 2, 1>{r, 0.5};
        if (data.n == 2 && data.p == 2 && data.q == 0) 
            return influence::exponential_2d<T, 2, 2>{r, 0.5};
        if (data.n == 2 && data.p == 3 && data.q == 0) 
            return influence::exponential_2d<T, 3, 2>{r, 0.5};
        if (data.n == 2 && data.p == 5 && data.q == 0) 
            return influence::exponential_2d<T, 5, 2>{r, 0.5};
        if (data.n == 5 && data.p == 2 && data.q == 0)
            return influence::exponential_2d<T, 2, 5>{r, 0.5};
        if (data.n == std::numeric_limits<size_t>::max() && data.p == 2 && data.q == 0)
            return influence::exponential_2d<T, 2, std::numeric_limits<size_t>::max()>{r, 0.5};

        if (data.n == 2 && data.p == 2 && data.q == 1)
            return influence::exponential_2d<T, 2>{r, T(data.q)};
        if (data.n == 2 && data.p == 2 && data.q == 2) // q == 2 temporary condition. In future refactoring will be changed 
            return influence::exponential_2d<T, 2>{r, 1.5};
        if (data.n == 2 && data.p == 2 && data.q == 3)
            return influence::exponential_2d<T, 2>{r, T(data.q)};

        if (data.n == 2 && data.p == 1 && data.q == 1)
            return influence::exponential_2d<T, 1>{r, T(data.q)};
        if (data.n == 2 && data.p == 3 && data.q == 1)
            return influence::exponential_2d<T, 3>{r, T(data.q)};

        return influence::exponential_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
    }
    throw std::domain_error{"Unsupported influence family."};
}

template<std::floating_point T, std::signed_integral I>
void save_csv(const std::optional<thermal::heat_equation_solution_2d<T, I>>& thermal_solution,
              const std::optional<mechanical::mechanical_solution_2d<T, I>>& mechanical_solution,
              const config::save_data& save) {
    if (parallel::MPI_rank() != 0 || !save.contains("csv")) // Only the master process saves data
        return;
    std::vector<std::pair<std::string, const std::vector<T>&>> data;
    if (thermal_solution) {
        data.push_back({"temperature", thermal_solution->temperature()});
        if (thermal_solution->is_flux_calculated()) {
            data.push_back({"flux_x", thermal_solution->flux()[X]});
            data.push_back({"flux_y", thermal_solution->flux()[Y]});
        }
    }
    if (mechanical_solution) {
        data.push_back({"displacement_x", mechanical_solution->displacement()[X]});
        data.push_back({"displacement_y", mechanical_solution->displacement()[Y]});
        if (mechanical_solution->is_strain_and_stress_calculated()) {
            data.push_back({"strain_11", mechanical_solution->strain()[0]});
            data.push_back({"strain_22", mechanical_solution->strain()[1]});
            data.push_back({"strain_12", mechanical_solution->strain()[2]});
            data.push_back({"stress_11", mechanical_solution->stress()[0]});
            data.push_back({"stress_22", mechanical_solution->stress()[1]});
            data.push_back({"stress_12", mechanical_solution->stress()[2]});
        }
    }
    if (data.empty())
        throw std::logic_error{"Nothig to save."};
    const auto& container = thermal_solution ? thermal_solution->mesh().container() : mechanical_solution->mesh().container();
    mesh::utils::save_as_csv(save.path("csv", "csv"), container, data, save.precision());
}

template<std::floating_point T, std::signed_integral I>
void save_vtk(const std::optional<thermal::heat_equation_solution_2d<T, I>>& thermal_solution,
              const std::optional<mechanical::mechanical_solution_2d<T, I>>& mechanical_solution,
              const config::save_data& save) {
    if (parallel::MPI_rank() != 0 || !save.contains("vtk")) // Only the master process saves data
        return;
    if (!thermal_solution && !mechanical_solution)
        throw std::logic_error{"Nothig to save."};
    std::ofstream vtk{save.path("vtk", "vtk")};
    vtk.precision(save.precision() ? *save.precision() : vtk.precision());
    const auto& container = thermal_solution ? thermal_solution->mesh().container() : mechanical_solution->mesh().container();
    mesh::utils::save_as_vtk(vtk, container);
    vtk << "POINT_DATA " << container.nodes_count() << '\n';
    if (thermal_solution) {
        mesh::utils::save_scalars_to_vtk(vtk, "temperature", thermal_solution->temperature());
        if (thermal_solution->is_flux_calculated())
            mesh::utils::save_vectors_to_vtk(vtk, "flux", thermal_solution->flux());
    }
    if (mechanical_solution) {
        mesh::utils::save_vectors_to_vtk(vtk, "displacement", mechanical_solution->displacement());
        if (mechanical_solution->is_strain_and_stress_calculated()) {
            mesh::utils::save_tensors_to_vtk(vtk, "strain", mechanical_solution->strain());
            mesh::utils::save_tensors_to_vtk(vtk, "stress", mechanical_solution->stress());
        }
    }
}

}
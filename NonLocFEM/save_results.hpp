#pragma once

#include <config/save_data.hpp>
#include <mesh/mesh_1d/mesh_1d.hpp>
#include <solvers/solver_2d/mechanical/mechanical_solution_2d.hpp>
#include <solvers/solver_2d/thermal/heat_equation_solution_2d.hpp>

namespace nonlocal {

template<std::floating_point T>
void save_csv(const std::optional<solver_2d::thermal::heat_equation_solution_2d<T>>& thermal_solution,
              const std::optional<solver_2d::mechanical::mechanical_solution_2d<T>>& mechanical_solution,
              const config::save_data& save, const std::optional<uint64_t> step = std::nullopt) {
    if (parallel::MPI_rank() != 0 || !save.contains("csv")) // Only the master process saves data
        return;
    using namespace std::string_literals;
    std::vector<mesh::utils::data_pair_t<T>> data;
    if (thermal_solution) {
        data.emplace_back(std::pair{"temperature"s, std::cref(thermal_solution->temperature())});
        if (thermal_solution->is_flux_calculated())
           data.emplace_back(std::pair{std::array{"flux_x"s, "flux_y"s}, std::cref(thermal_solution->flux())});
    }
    if (mechanical_solution) {
        data.emplace_back(std::pair{std::array{"displacement_x"s, "displacement_y"s}, std::cref(mechanical_solution->displacement())});
        if (mechanical_solution->is_strain_and_stress_calculated()) {
            data.emplace_back(std::pair{std::array{"strain_11"s, "strain_22"s, "strain_12"s}, std::cref(mechanical_solution->strain())});
            data.emplace_back(std::pair{std::array{"stress_11"s, "stress_22"s, "stress_12"s}, std::cref(mechanical_solution->stress())});
        }
    }
    if (data.empty())
        throw std::logic_error{"Nothig to save."};
    const auto& container = thermal_solution ? thermal_solution->mesh().container() : mechanical_solution->mesh().container();
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + "_" + save.get_name("csv"), "csv") : 
                                              save.path("csv", "csv");
    mesh::utils::save_as_csv(path, container, data, save.precision());
}

template<std::floating_point T>
void save_vtk(const std::optional<solver_2d::thermal::heat_equation_solution_2d<T>>& thermal_solution,
              const std::optional<solver_2d::mechanical::mechanical_solution_2d<T>>& mechanical_solution,
              const config::save_data& save, const std::optional<uint64_t> step = std::nullopt) {
    if (parallel::MPI_rank() != 0 || !save.contains("vtk")) // Only the master process saves data
        return;
    if (!thermal_solution && !mechanical_solution)
        throw std::logic_error{"Nothig to save."};
    const std::filesystem::path path = step ? save.make_path(std::to_string(*step) + "_" + save.get_name("vtk"), "vtk") : 
                                              save.path("vtk", "vtk", "solution");
    std::ofstream vtk{path};
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
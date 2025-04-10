#pragma once

#include <config/save_data.hpp>
#include <mesh/mesh_1d/mesh_1d.hpp>
#include <solvers/solver_2d/mechanical/mechanical_solution_2d.hpp>
#include <solvers/solver_2d/thermal/heat_equation_solution_2d.hpp>

namespace nonlocal {

template<std::floating_point T, std::integral I>
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

template<std::floating_point T, std::integral I>
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
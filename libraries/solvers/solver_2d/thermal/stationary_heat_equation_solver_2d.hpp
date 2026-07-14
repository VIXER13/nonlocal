#pragma once

#include "conductivity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "radiation_condition_2d.hpp"
#include "init_problem_settings.hpp"
#include "heat_equation_solution_2d.hpp"

#include <solvers/base/utils.hpp>
#include <solvers/slae/conjugate_gradient.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

namespace nonlocal::solver_2d::thermal {

template<class T>
struct stationary_equation_parameters_2d final {
    std::optional<std::function<T(const std::array<T, 2>&)>> right_part;
    std::optional<std::function<T(const std::array<T, 2>&)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-8;
    size_t max_iterations = 100;
    T energy = T{0};
};

template<class T, class I>
std::vector<T> init_right_part(const mesh::mesh_2d<T, I>& mesh,
                               const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                               const stationary_equation_parameters_2d<T>& auxiliary_data,
                               const bool is_neumann = false) {
    static constexpr size_t DoF = 1;
    std::vector<T> right_part(mesh.container().nodes_count() + is_neumann, T{0});
    boundary_condition_second_kind_2d(right_part, mesh, boundaries_conditions);
    if (auxiliary_data.right_part)
        integrate_right_part<DoF>(right_part, mesh, *auxiliary_data.right_part);                                 
    if (is_neumann) {
        if (std::abs(std::reduce(right_part.begin(), right_part.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
            throw std::domain_error{"It's unsolvable Neumann problem."};
        right_part[right_part.size() - 1] = auxiliary_data.energy;
    }
    return right_part;
}

template<class Matrix_Index, class T, class I>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                   const parameters_2d<T>& parameters,
                                                                   const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                   const stationary_equation_parameters_2d<T>& auxiliary_data) {
    using namespace metamath::operators;
    static constexpr bool Is_Stationary = true;
    const auto settings = init_problem_settings(mesh->container(), parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);

    std::vector<T> f_init = init_right_part(*mesh, boundaries_conditions, auxiliary_data, settings.is_neumann); 
    std::vector<T> f = f_init;
    std::vector<T> temperature_curr(mesh->container().nodes_count() + settings.is_neumann, T{0});
    if (auxiliary_data.initial_distribution)
        for(const size_t node : mesh->container().nodes())
            temperature_curr[node] = (*auxiliary_data.initial_distribution)(mesh->container().node_coord(node));
    std::vector<T> temperature_prev = temperature_curr;
    std::vector<T> residual(mesh->container().nodes_count() + settings.is_neumann, T{0});
    std::vector<T> residual_rad(mesh->container().nodes_count() + settings.is_neumann, T{0});

    const auto conductivity_parameters = evaluate_conductivity(*mesh, parameters, std::vector<T>(mesh->quad_shift(mesh->container().elements_2d_count()), T{0}));
    T difference = T{1};
    T norm_of_residual = T{1};
    size_t iteration = 0;
    while (iteration < auxiliary_data.max_iterations && 
           auxiliary_data.tolerance < difference && 
           auxiliary_data.tolerance < norm_of_residual) {
        std::swap(temperature_curr, temperature_prev);
        f = f_init;
        conductivity_matrix_2d<T, I> conductivity{mesh};
        conductivity.compute(conductivity_parameters, settings.is_inner_nodes, settings.is_symmetric(), settings.is_neumann);
        convection_condition_2d(conductivity.matrix().inner(), *mesh, boundaries_conditions);
        residual = conductivity.matrix().inner().template self_adjoint<metamath::linear::matrix_part::Upper>() * temperature_prev;
        std::fill(residual_rad.begin(), residual_rad.end(), T{0});
        radiation_condition_2d(conductivity.matrix().inner(), residual_rad, *mesh, boundaries_conditions, temperature_prev);
        if (!settings.is_neumann) {
            first_kind_matrix_fill_2d(conductivity.matrix().inner(), residual_rad, *mesh, boundaries_conditions);
            boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, conductivity.matrix().bound());
        }           
        residual -= residual_rad + f;
        if (settings.is_symmetric()) {
            logger::info() << "Local matrix" << std::endl;
            if (settings.is_nonlinear_boundary) {                   
                //const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{conductivity.matrix().inner()};
                //temperature_curr = temperature_prev - solver.solve(residual);
            } else {
                // conductivity_matrix_2d<T, I> conductivity_local{mesh};
                // conductivity_local.nodes_for_processing(std::ranges::iota_view<size_t, size_t>{0u, mesh->container().nodes_count()});
                // conductivity_local.compute(conductivity_parameters, settings.is_inner_nodes, settings.is_symmetric(), settings.is_neumann, assemble_part::LOCAL);
                // slae::conjugate_gradient<T> solver{conductivity.matrix().inner()};
                // logger::info() << "ILLT preconditioner" << std::endl;
                // solver.template init_preconditioner<slae::eigen_ILLT_preconditioner>(conductivity_local.matrix().inner());
                // if (solver.preconditioner().computation_info() != Eigen::Success) {
                //     solver.template init_preconditioner<slae::eigen_identity_preconditioner>();
                //     logger::warning() << "The ILLT preconditioner could not be calculated, "
                //                       << "the preconditioner was switched to Identity." << std::endl;
                // } 
                // temperature_curr = temperature_prev - solver.solve(residual);  
            }
        } else {
            //const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>> solver{conductivity.matrix().inner()};
            //temperature_curr = temperature_prev - solver.solve(residual);
        }
        if (!settings.is_nonlinear_boundary)
            break;
        ++iteration;
        norm_of_residual = metamath::linear::norm(residual);
        difference = metamath::linear::norm(temperature_curr - temperature_prev) / (metamath::linear::norm(temperature_curr) ?: T{1});
        logger::info() << " Iteration #"           << iteration 
                       << ", norm(prev - curr) = " << difference 
                       << ", residual = "          << norm_of_residual << ";" << std::endl;
    }
    heat_equation_solution_2d<T, I> solution{mesh, conductivity_parameters, temperature_curr};
    solution.calc_flux();
    return solution;
}

}
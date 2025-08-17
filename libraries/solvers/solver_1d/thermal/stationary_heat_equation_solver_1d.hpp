#pragma once

#include "thermal_conductivity_assembler.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "init_problem_settings.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_1d final {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T energy = T{0};
};

template<class T>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_right_part(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                    const thermal_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                    const stationary_equation_parameters_1d<T>& additional_parameters,
                                                    const bool is_neumann = false) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(right_part, boundaries_conditions, is_neumann);
    if (additional_parameters.right_part)
        integrate_right_part(right_part, *mesh, *additional_parameters.right_part);
    if (is_neumann && std::abs(std::reduce(right_part.begin(), right_part.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};
    if (is_neumann)
        right_part[right_part.size() - 1] = additional_parameters.energy;
    return right_part;
}

template<class T, class I>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const parameters_1d<T>& parameters,
                                                                const thermal_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                                const stationary_equation_parameters_1d<T>& additional_parameters) {
    static constexpr bool Stationary_Problem = true;
    const auto settings = init_problem_settings(parameters, boundaries_conditions, Stationary_Problem);
    log_problem_settings(settings);

    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = init_right_part(mesh, boundaries_conditions, additional_parameters, settings.is_neumann);
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_right_part = right_part;
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(right_part.size());

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(right_part.size());
    if (settings.is_nonlinear() && additional_parameters.initial_distribution)
        for(const size_t node : mesh->nodes())
            temperature_prev[node] = (*additional_parameters.initial_distribution)(mesh->node_coord(node));
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = temperature_prev;
    
    finite_element_matrix_1d<T, I> conductivity;
    init_matrix_portrait(conductivity.inner, *mesh, settings);
    thermal_conductivity_assembler_1d<T, I> assembler{conductivity, mesh};

    T difference = T{1};
    size_t iteration = 0;
    do {
        if (settings.is_nonlinear()) {
            std::swap(temperature_prev, temperature_curr);
            std::copy(initial_right_part.begin(), initial_right_part.end(), right_part.begin());
            conductivity.set_zero();
        }
        using nonlocal::mesh::utils::from_nodes_to_qnodes;
        assembler.calc_matrix(parameters, settings,
            settings.is_solution_dependent ? std::optional{from_nodes_to_qnodes(*mesh, temperature_prev)} : std::nullopt
        );

        if (settings.is_neumann) {
            if (settings.is_symmetric()) {
                const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.inner};
                temperature_curr = solver.solveWithGuess(right_part, temperature_prev);
            } else {
                const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, I>> solver{conductivity.inner};
                temperature_curr = solver.solveWithGuess(right_part, temperature_prev);
            }
        } else {
            convection_condition_1d(conductivity.inner, boundaries_conditions);
            boundary_condition_first_kind_1d(right_part, conductivity.bound, boundaries_conditions);
            if (settings.is_radiation_boundary) {
                if (settings.is_symmetric())
                    residual = conductivity.inner.template selfadjointView<Eigen::Upper>() * temperature_prev - right_part;
                else
                    residual = conductivity.inner * temperature_prev - right_part;
                radiation_condition_1d(conductivity.inner, residual, boundaries_conditions, temperature_prev);
            }

            if (settings.is_symmetric()) {
                const Eigen::SimplicialCholesky<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>
                > solver{conductivity.inner};
                if (settings.is_radiation_boundary)
                    temperature_curr = temperature_prev - solver.solve(residual);
                else
                    temperature_curr = solver.solve(right_part);
            } else {
                const Eigen::SparseLU<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::NaturalOrdering<I>
                > solver{conductivity.inner};
                if (settings.is_radiation_boundary)
                    temperature_curr = temperature_prev - solver.solve(residual);
                else
                    temperature_curr = solver.solve(right_part);
            }
        }
        ++iteration;
        difference = (temperature_curr - temperature_prev).norm() / (temperature_curr.norm() ?: T{1});
    } while(settings.is_nonlinear() &&
            iteration < additional_parameters.max_iterations && 
            difference > additional_parameters.tolerance);
    return heat_equation_solution_1d<T>{mesh, parameters, temperature_curr};
}

}
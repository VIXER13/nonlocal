#pragma once

#include "conductivity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "radiation_condition_2d.hpp"
#include "thermal_boundary_conditions_2d.hpp"
#include "heat_equation_solution_2d.hpp"

#include <solvers/solvers_utils.hpp>
#include <solvers/slae/conjugate_gradient.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

#include <chrono>

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_2d final {
    std::optional<std::function<T(const std::array<T, 2>&)>> right_part;
    std::optional<std::function<T(const std::array<T, 2>&)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-8;
    size_t max_iterations = 100;
    T energy = T{0};
};

template<class T>
bool is_neumann_problem(const thermal_boundaries_conditions_2d<T>& boundaries_conditions) {
    const auto conditions = boundaries_conditions | std::views::values;
    return std::all_of(conditions.begin(), conditions.end(), [](const auto& condition) {
        return bool(dynamic_cast<const flux_2d<T>*>(condition.get()));
    });
}

template<class T>
bool is_radiation_problem(const thermal_boundaries_conditions_2d<T>& boundaries_conditions) {
    const auto conditions = boundaries_conditions | std::views::values;
    return std::any_of(conditions.begin(), conditions.end(), [](const auto& condition) {
        return bool(dynamic_cast<const radiation_2d<T>*>(condition.get()));
    });
}

template<class T>
bool is_nonlinear_problem(const parameters_2d<T>& parameters) {
    return false;
    // return std::any_of(parameters.begin(), parameters.end(), [](const auto& parameter) {
    //     return parameter.second.physical->type != coefficients_t::CONSTANTS;
    // });
}

template<class T>
bool is_nonlocal_problem(const parameters_2d<T>& parameters) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    return std::any_of(theories.begin(), theories.end(), [](const auto& theory) noexcept {
        return theory.second == theory_t::NONLOCAL;
    });
}

template<class Matrix_Index, class T, class I>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                   const parameters_2d<T>& parameters,
                                                                   const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                   const stationary_equation_parameters_2d<T>& auxiliary_data) {
    static constexpr size_t DoF = 1;
    const bool is_neumann = is_neumann_problem(boundaries_conditions);
    const bool is_radiation = is_radiation_problem(boundaries_conditions);
    if (is_neumann)
        logger::info() << "Neumann problem" << std::endl;
    const bool is_symmetric = !(is_nonlinear_problem(parameters) && is_nonlocal_problem(parameters));
    logger::info() << (is_symmetric ? "Symmetric" : "Asymmetrical") << " problem" << std::endl;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count() + is_neumann);
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    if (auxiliary_data.right_part)
        integrate_right_part<DoF>(f, *mesh, *auxiliary_data.right_part);
    if (is_neumann) {
    //     if (!is_solvable_neumann_problem(*mesh_proxy, f))
    //         throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
        f[f.size() - 1] = auxiliary_data.energy;
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count() + is_neumann);
    if (auxiliary_data.initial_distribution)
        for(const size_t node : mesh->container().nodes())
            temperature_curr[node] = (*auxiliary_data.initial_distribution)(mesh->container().node_coord(node));
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = temperature_curr;
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual = temperature_curr;

    T difference = T{1}, norm_of_res = T{1};
    size_t iteration = 0;
    while (iteration < auxiliary_data.max_iterations && 
           auxiliary_data.tolerance < difference && 
           auxiliary_data.tolerance < norm_of_res) {
        std::swap(temperature_curr, temperature_prev);
        conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh};
        conductivity.compute(parameters, utils::inner_nodes(mesh->container(), boundaries_conditions), is_symmetric, is_neumann);
        convection_condition_2d(conductivity.matrix().inner(), *mesh, boundaries_conditions);
        if (!is_neumann)
            boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, conductivity.matrix().bound());
        residual = conductivity.matrix().inner().template selfadjointView<Eigen::Upper>() * temperature_prev - f;
        radiation_condition_2d(conductivity.matrix().inner(), residual, *mesh, boundaries_conditions, temperature_prev);
        if (is_symmetric) {
            logger::info() << "Local matrix" << std::endl;
            conductivity_matrix_2d<T, I, Matrix_Index> conductivity_local{mesh};
            conductivity_local.nodes_for_processing(std::ranges::iota_view<size_t, size_t>{0u, mesh->container().nodes_count()});
            conductivity_local.compute(parameters, utils::inner_nodes(mesh->container(), boundaries_conditions), is_symmetric, is_neumann, assemble_part::LOCAL);
            slae::conjugate_gradient<T, Matrix_Index> solver{conductivity.matrix().inner()};
            logger::info() << "ILLT preconditioner" << std::endl;
            solver.template init_preconditioner<slae::eigen_ILLT_preconditioner>(
                conductivity_local.matrix().inner()
            );
            if (solver.preconditioner().computation_info() != Eigen::Success) {
                solver.template init_preconditioner<slae::eigen_identity_preconditioner>();
                logger::warning() << "The ILLT preconditioner could not be calculated, "
                                << "the preconditioner was switched to Identity." << std::endl;
            }
            temperature_curr = temperature_prev - solver.solve(residual);
        } else {
            const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>> solver{conductivity.matrix().inner()};
            temperature_curr = temperature_prev - solver.solve(residual);
        }

        if (!is_radiation)
            break;
        ++iteration;
        norm_of_res = residual.norm();
        difference = (temperature_curr - temperature_prev).norm() / (temperature_curr.norm() ?: T{1});
        logger::info() << "Iteration #" << iteration <<", norm(prev - curr) = " << difference << ", residual = " << norm_of_res << ";" << std::endl;
    }
    return heat_equation_solution_2d<T, I>{mesh, parameters, temperature_curr};
}

}
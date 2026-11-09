#pragma once

#include "assemble_matrix_portrait.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "linear_solver_1d.hpp"
#include "right_part_1d.hpp"
#include "spatial_assembler_1d.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>

#include <numeric>

namespace nonlocal::solver_1d {

template<std::floating_point T>
struct stationary_solver_parameters_1d final {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T integral_value = T{ 0 }; // Pure Neumann problem: integral(u(x) dx) = integral_value.
};

template<std::floating_point T, physics_t Physics>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_stationary_right_part_1d(
    const mesh::mesh_1d<T>& mesh, const boundaries_conditions_1d<T, Physics>& boundaries,
    const stationary_solver_parameters_1d<T>& parameters, const bool is_neumann = false) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh.nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(right_part, boundaries, is_neumann);
    if (parameters.right_part)
        integrate_right_part(right_part, mesh, *parameters.right_part);
    if (is_neumann) {
        if (std::abs(std::reduce(right_part.begin(), right_part.end(), T{ 0 })) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
            throw std::domain_error{ "It's unsolvable Neumann problem." };
        right_part[right_part.size() - 1] = parameters.integral_value;
    }
    return right_part;
}

template<std::floating_point T, std::integral I, physics_t Physics, class Parameters, class Get_Coefficient>
Eigen::Matrix<T, Eigen::Dynamic, 1> stationary_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const Parameters& parameters,
    const boundaries_conditions_1d<T, Physics>& boundaries, const problem_settings& settings,
    const stationary_solver_parameters_1d<T>& solver_parameters, const Get_Coefficient get_coefficient,
    void (*apply_robin)(Eigen::SparseMatrix<T, Eigen::RowMajor, I>&,
                        const boundaries_conditions_1d<T, Physics>&) = nullptr,
    void (*apply_nonlinear)(Eigen::SparseMatrix<T, Eigen::RowMajor, I>&, Eigen::Matrix<T, Eigen::Dynamic, 1>&,
                            const boundaries_conditions_1d<T, Physics>&, const Eigen::Matrix<T, Eigen::Dynamic, 1>&,
                            T) = nullptr) {
    if (settings.is_nonlinear_boundary && !apply_nonlinear)
        throw std::invalid_argument{ "A nonlinear boundary handler is required." };
    log_problem_settings(settings);
    using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    const vector initial_right_part = init_stationary_right_part_1d(*mesh, boundaries, solver_parameters, settings.is_neumann);
    vector solution_prev = vector::Zero(initial_right_part.size());
    if (settings.is_nonlinear() && solver_parameters.initial_distribution)
        for (const size_t node : mesh->nodes())
            solution_prev[node] = (*solver_parameters.initial_distribution)(mesh->node_coord(node));
    vector solution_curr = solution_prev;

    finite_element_matrix_1d<T, I> matrix;
    init_matrix_portrait(matrix.inner, *mesh, settings);
    spatial_assembler_1d<T, I> assembler{ matrix, mesh };
    // The Neumann multiplier and a nonlinear boundary Jacobian can make a symmetric system indefinite.
    const linear_system_properties_1d properties{
        .is_symmetric = settings.is_symmetric(),
        .may_be_indefinite = settings.is_neumann || settings.is_nonlinear_boundary
    };

    T difference = T{ 1 };
    size_t iteration = 0;
    do {
        if (settings.is_nonlinear()) {
            std::swap(solution_prev, solution_curr);
            matrix.set_zero();
        }
        vector right_part = initial_right_part;
        using nonlocal::mesh::utils::from_nodes_to_qnodes;
        assembler.calc_matrix(parameters, settings, get_coefficient,
                              settings.is_solution_dependent ? std::optional{ from_nodes_to_qnodes(*mesh, solution_prev) }
                                                             : std::nullopt);
        matrix.set_first_kind_diagonal(settings.is_first_kind);
        if (!settings.is_neumann) {
            if (apply_robin)
                apply_robin(matrix.inner, boundaries);
            boundary_condition_first_kind_1d(right_part, matrix.bound, boundaries);
        }

        if (settings.is_nonlinear_boundary) {
            vector residual(right_part.size());
            if (properties.is_symmetric)
                residual = matrix.inner.template selfadjointView<Eigen::Upper>() * solution_prev - right_part;
            else
                residual = matrix.inner * solution_prev - right_part;
            apply_nonlinear(matrix.inner, residual, boundaries, solution_prev, T{ 1 });
            solution_curr = solution_prev - solve_linear_system_1d(matrix.inner, residual, properties);
        } else
            solution_curr = solve_linear_system_1d(matrix.inner, right_part, properties);

        ++iteration;
        const T norm = solution_curr.norm();
        difference = (solution_curr - solution_prev).norm() / (norm == T{ 0 } ? T{ 1 } : norm);
    } while (settings.is_nonlinear() && iteration < solver_parameters.max_iterations && difference > solver_parameters.tolerance);
    return solution_curr;
}

} // namespace nonlocal::solver_1d

#pragma once

#include "init_problem_settings.hpp"
#include "mass_matrix_assembler.hpp"
#include "mechanical_equation_solution_1d.hpp"
#include "spring_condition_1d.hpp"
#include "stiffness_matrix_assembler.hpp"

#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/linear_solver_1d.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
struct time_harmonic_equation_parameters_1d {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T frequency = T{ 0 };
};

template<std::floating_point T>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_right_part(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                    const mechanical_boundaries_conditions_1d<T>& boundaries_conditions,
                                                    const time_harmonic_equation_parameters_1d<T>& additional_parameters,
                                                    const bool is_neumann = false) {
    return assemble_right_part_1d(*mesh, boundaries_conditions, additional_parameters.right_part, is_neumann);
}

template<std::floating_point T, std::integral I>
mechanical_equation_solution_1d<T> harmonic_mechanical_equation_solver_1d(
    const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const parameters_1d<T>& parameters,
    const mechanical_boundaries_conditions_1d<T>& boundaries_conditions,
    const time_harmonic_equation_parameters_1d<T>& additional_parameters = {}) {
    static constexpr bool Is_Stationary = false;
    const auto settings = init_problem_settings(parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);

    finite_element_matrix_1d<T, I> stiffness;
    init_matrix_portrait(stiffness.inner, *mesh, settings);
    stiffness_assembler_1d<T, I> stiffness_assembler{ stiffness, mesh };
    stiffness_assembler.calc_matrix(parameters, settings);

    const problem_settings mass_settings{
        .theories = std::vector<theory_t>(mesh->segments_count(), theory_t::LOCAL),
        .is_first_kind = settings.is_first_kind
    };
    finite_element_matrix_1d<T, I> mass;
    init_matrix_portrait(mass.inner, *mesh, mass_settings);
    mass_assembler_1d<T, I> mass_assembler{ mass, mesh };
    mass_assembler.calc_matrix(parameters, settings.is_first_kind);

    const T omega_square = additional_parameters.frequency * additional_parameters.frequency;
    const bool is_symmetric = settings.is_symmetric();
    if (is_symmetric)
        stiffness.inner -= omega_square * mass.inner;
    else {
        // Mass is stored as an upper triangle. Subtract both halves when stiffness is stored in full.
        const Eigen::SparseMatrix<T, Eigen::RowMajor, I> full_mass = mass.inner.template selfadjointView<Eigen::Upper>();
        stiffness.inner -= omega_square * full_mass;
    }
    for (const size_t b : std::ranges::iota_view{ 0u, 2u })
        for (const auto& [col, value] : mass.bound[b])
            stiffness.bound[b][col] -= omega_square * value;

    stiffness.set_first_kind_diagonal(settings.is_first_kind);

    auto right_part = init_right_part(mesh, boundaries_conditions, additional_parameters);
    spring_condition_1d(stiffness.inner, boundaries_conditions);
    boundary_condition_first_kind_1d(right_part, stiffness.bound, boundaries_conditions);

    // K - omega^2*M can be indefinite even when it is symmetric.
    const auto displacement = solve_linear_system_1d(stiffness.inner, right_part,
        { .is_symmetric = is_symmetric, .may_be_indefinite = omega_square != T{ 0 } });
    return mechanical_equation_solution_1d<T>{ mesh, parameters, displacement };
}

} // namespace nonlocal::solver_1d::mechanical

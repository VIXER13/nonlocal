#pragma once

#include "stiffness_matrix_assembler.hpp"
#include "mass_matrix_assembler.hpp"
#include "mechanical_equation_solution_1d.hpp"
#include "spring_condition_1d.hpp"
#include "init_problem_settings.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
struct time_harmonic_equation_parameters_1d {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T frequency = T{0};
};

template<std::floating_point T>
Eigen::Matrix<T, Eigen::Dynamic, 1> init_right_part(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                    const mechanical_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                    const time_harmonic_equation_parameters_1d<T>& additional_parameters,
                                                    const bool is_neumann = false) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(right_part, boundaries_conditions, is_neumann);
    if (additional_parameters.right_part)
        integrate_right_part(right_part, *mesh, *additional_parameters.right_part);
    if (is_neumann && std::abs(std::reduce(right_part.begin(), right_part.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};
    return right_part;
}

template<std::floating_point T, std::integral I>
mechanical_equation_solution_1d<T> harmonic_mechanical_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                          const parameters_1d<T>& parameters,
                                                                          const mechanical_boundaries_conditions_1d<T>& boundaries_conditions,
                                                                          const time_harmonic_equation_parameters_1d<T>& additional_parameters = {}) {
    static constexpr bool Is_Stationary = false; 
    const auto settings = init_problem_settings(parameters, boundaries_conditions, Is_Stationary);
    log_problem_settings(settings);

    finite_element_matrix_1d<T, I> stiffness;
    init_matrix_portrait(stiffness.inner, *mesh, settings);
    stiffness_assembler_1d<T, I> stiffness_assembler{stiffness, mesh};
    stiffness_assembler.calc_matrix(parameters, settings);

    finite_element_matrix_1d<T, I> mass;
    init_matrix_portrait(mass.inner, *mesh, settings);
    mass_assembler_1d<T, I> mass_assembler{mass, mesh};
    mass_assembler.calc_matrix(parameters, settings.is_first_kind);

    const T omega_square = additional_parameters.frequency * additional_parameters.frequency;
    stiffness.inner -= omega_square * mass.inner;
    for(const size_t b : std::ranges::iota_view{0u, 2u})
        for(const auto& [col, value] : mass.bound[b])
            stiffness.bound[b][col] -= omega_square * value;

    // The first kind rows of both K and M carry a unit diagonal; after the
    // subtraction the Dirichlet diagonal becomes (1 - w^2). Reset it to 1.
    const size_t last_node = mesh->nodes_count() - 1;
    if (settings.is_first_kind.front())
        stiffness.inner.coeffRef(0, 0) = T{1};
    if (settings.is_first_kind.back())
        stiffness.inner.coeffRef(last_node, last_node) = T{1};

    Eigen::Matrix<T, Eigen::Dynamic, 1> right_part = init_right_part(mesh, boundaries_conditions, additional_parameters, false);
    spring_condition_1d(stiffness.inner, boundaries_conditions);
    boundary_condition_first_kind_1d(right_part, stiffness.bound, boundaries_conditions);

    Eigen::Matrix<T, Eigen::Dynamic, 1> displacement;
    if (settings.is_symmetric()) {
        const Eigen::SimplicialCholesky<
            Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>
        > solver{stiffness.inner};
        displacement = solver.solve(right_part);
    } else {
        const Eigen::SparseLU<
            Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::NaturalOrdering<I>
        > solver{stiffness.inner};
        displacement = solver.solve(right_part);
    }
    return mechanical_equation_solution_1d<T>{mesh, parameters, displacement};
}

}

#pragma once

#include "solvers_utils.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "thermal_boundary_conditions_2d.hpp"
#include "boundary_condition_first_kind_2d.hpp"
#include "boundary_condition_second_kind_2d.hpp"
#include "convection_condition_2d.hpp"
#include "right_part_2d.hpp"
#include "heat_equation_solution_2d.hpp"

#include "conjugate_gradient.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_2d final {
    std::optional<std::function<T(const std::array<T, 2>&)>> right_part;
    std::optional<std::function<T(const std::array<T, 2>&)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
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
bool is_nonlinear_problem(const parameters_2d<T>& parameters) {
    return std::any_of(parameters.begin(), parameters.end(), [](const auto& parameter) {
        return parameter.second.physical->type != coefficients_t::CONSTANTS;
    });
}

template<class T>
bool is_nonlocal_problem(const parameters_2d<T>& parameters) {
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    return std::any_of(theories.begin(), theories.end(), [](const auto& theory) noexcept {
        return theory.second == theory_t::NONLOCAL;
    });
}

template<class Matrix_Index, class T, class I, class Right_Part>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                   const parameters_2d<T>& parameters,
                                                                   const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                   const Right_Part& right_part,
                                                                   const T energy = T{0}) {
    static constexpr size_t DoF = 1;
    const bool is_neumann = is_neumann_problem(boundaries_conditions);
    if (is_neumann)
        logger::info() << "Neumann problem" << std::endl;
    const bool is_symmetric = !(is_nonlinear_problem(parameters) && is_nonlocal_problem(parameters));
    logger::info() << (is_symmetric ? "Symmetric" : "Asymmetrical") << " problem" << std::endl;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count() + is_neumann);
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    integrate_right_part<DoF>(f, *mesh, right_part);
    if (is_neumann) {
    //     if (!is_solvable_neumann_problem(*mesh_proxy, f))
    //         throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
        f[f.size() - 1] = energy;
    }

    logger::info() << "Original matrix" << std::endl;
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh};
    conductivity.compute(parameters, utils::inner_nodes(mesh->container(), boundaries_conditions), is_symmetric, is_neumann);
    convection_condition_2d(conductivity.matrix().inner(), *mesh, boundaries_conditions);
    if (!is_neumann)
       boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, conductivity.matrix().bound());

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature;
    if (is_symmetric) {
        logger::info() << "Local matrix" << std::endl;
        thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity_local{mesh};
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
        temperature = solver.solve(f);
    } else {
        const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>> solver{conductivity.matrix().inner()};
        temperature = solver.solve(f);
    }
    return heat_equation_solution_2d<T, I>{mesh, parameters, temperature};
}

}
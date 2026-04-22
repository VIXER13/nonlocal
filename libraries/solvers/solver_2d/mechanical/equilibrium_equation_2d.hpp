#pragma once

#include "evaluate_mechanical_parameters.hpp"
#include "stiffness_matrix_2d.hpp"
#include "mechanical_boundary_conditions_2d.hpp"
#include "mechanical_solution_2d.hpp"
#include "init_problem_settings.hpp"
#include "temperature_condition_2d.hpp"

#include <solvers/base/utils.hpp>
#include <solvers/slae/conjugate_gradient.hpp>
#include <solvers/solver_2d/base/boundary_condition_first_kind_2d.hpp>
#include <solvers/solver_2d/base/boundary_condition_second_kind_2d.hpp>
#include <solvers/solver_2d/base/right_part_2d.hpp>

#include <optional>

namespace nonlocal::solver_2d::mechanical {

template<class Matrix_Index, class T, class I>
mechanical::mechanical_solution_2d<T, I> equilibrium_equation(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                              const raw_mechanical_parameters<T>& parameters,
                                                              const mechanical_boundaries_conditions_2d<T>& boundaries_conditions,
                                                              const std::vector<T>& delta_temperature = {},
                                                              const std::optional<std::function<std::array<T, 2>(const std::array<T, 2>&)>>& right_part = std::nullopt) {
    const auto settings = init_problem_settings(mesh->container(), parameters, boundaries_conditions);
    log_problem_settings(settings);
    const auto evaluated_parameters = evaluate_mechanical_parameters(*mesh, parameters, delta_temperature);
    stiffness_matrix<T, I, Matrix_Index> stiffness{mesh};
    stiffness.compute(evaluated_parameters, settings);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(stiffness.matrix().inner().cols());
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    if (right_part)
        integrate_right_part<2>(f, *mesh, *right_part);
    temperature_condition(f, *mesh, evaluated_parameters);
    boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, stiffness.matrix().bound());
    stiffness_matrix<T, I, Matrix_Index> local_stiffness{mesh};
    if (settings.is_nonlocal()) {
        local_stiffness.nodes_for_processing(std::ranges::iota_view<size_t, size_t>{0u, mesh->container().nodes_count()});
        local_stiffness.compute(evaluated_parameters, settings, assemble_part::LOCAL);
    }

    Eigen::Matrix<T, Eigen::Dynamic, 1> displacement;
    if (settings.is_symmetric()) {
        slae::conjugate_gradient<T, Matrix_Index> solver{stiffness.matrix().inner()};
        if (settings.is_nonlocal()) {
            solver.template init_preconditioner<slae::eigen_ILLT_preconditioner>(local_stiffness.matrix().inner());
            if (solver.preconditioner().computation_info() != Eigen::Success) {
                solver.template init_preconditioner<slae::eigen_identity_preconditioner>();
                logger::warning() << "The ILLT preconditioner could not be calculated, "
                                << "the preconditioner was switched to Identity." << std::endl;
            }
        }
        displacement = solver.solve(f);
    } else {
        const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>> solver{stiffness.matrix().inner()};
        displacement = solver.solve(f);
    }
    auto solution = mechanical_solution_2d<T, I>{mesh, evaluated_parameters, displacement};
    solution.calc_strain_and_stress();
    return solution;
}

}
#pragma once

#include "solvers_utils.hpp"
#include "stiffness_matrix_2d.hpp"
#include "mechanical_boundary_conditions_2d.hpp"
#include "boundary_condition_first_kind_2d.hpp"
#include "boundary_condition_second_kind_2d.hpp"
#include "right_part_2d.hpp"
#include "mechanical_solution_2d.hpp"
#include "temperature_condition_2d.hpp"

#include "conjugate_gradient.hpp"

#include <chrono>

namespace nonlocal::mechanical {

template<class Matrix_Index, class T, class I, class Right_Part>
mechanical::mechanical_solution_2d<T, I> equilibrium_equation(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                              const mechanical_parameters_2d<T>& parameters,
                                                              const mechanical_boundaries_conditions_2d<T>& boundaries_conditions,
                                                              const Right_Part& right_part) {
    stiffness_matrix<T, I, Matrix_Index> stiffness{mesh};
    stiffness.compute(parameters.materials, parameters.plane, utils::inner_nodes(mesh->container(), boundaries_conditions));
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(stiffness.matrix().inner().cols());
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    integrate_right_part<2>(f, *mesh, right_part);
    temperature_condition(f, *mesh, parameters);
    stiffness_matrix<T, I, Matrix_Index> local_stiffness{mesh};
    local_stiffness.nodes_for_processing(std::ranges::iota_view<size_t, size_t>{0u, mesh->container().nodes_count()});
    local_stiffness.compute(parameters.materials, parameters.plane, utils::inner_nodes(mesh->container(), boundaries_conditions), assemble_part::LOCAL);
    slae::conjugate_gradient<T, Matrix_Index> solver{stiffness.matrix().inner()};
    solver.template init_preconditioner<slae::eigen_ILLT_preconditioner>(
        local_stiffness.matrix().inner()
    );
    if (solver.preconditioner().computation_info() != Eigen::Success) {
        solver.template init_preconditioner<slae::eigen_identity_preconditioner>();
        logger::get().log(logger::log_level::WARNING) << "The ILLT preconditioner could not be calculated, "
                                                      << "the preconditioner was switched to Identity." << std::endl;
    }
    const auto displacement = solver.solve(f);
    return mechanical_solution_2d<T, I>{mesh, parameters, displacement};
}

}
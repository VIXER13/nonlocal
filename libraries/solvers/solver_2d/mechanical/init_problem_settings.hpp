#pragma once

#include "mechanical_parameters_2d.hpp"
#include "mechanical_boundary_conditions_2d.hpp"

#include <solvers/solver_2d/base/problem_settings.hpp>
#include <solvers/solver_2d/base/solvers_utils.hpp>

#include <concepts>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T, std::integral I>
problem_settings init_problem_settings(const mesh::mesh_container_2d<T, I>& mesh,
                                       const raw_mechanical_parameters<T>& parameters,
                                       const mechanical_boundaries_conditions_2d<T>& boundaries_conditions) {
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) { return !is_constant(parameter.physical.elastic); };
    const auto parameters_view = parameters | std::views::values;
    return {
        .theories = theories_types(parameters),
        .is_nonconstant_parameters = std::any_of(parameters_view.begin(), parameters_view.end(), is_nonconstant_parameters),
        .is_inner_nodes = utils::inner_nodes(mesh, boundaries_conditions)
    };
}

}
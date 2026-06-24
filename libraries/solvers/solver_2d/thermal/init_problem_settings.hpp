#pragma once

#include "thermal_boundary_conditions_2d.hpp"
#include "thermal_parameters_2d.hpp"

#include <solvers/solver_2d/base/problem_settings.hpp>
#include <metamath/types/visitor.hpp>

#include <concepts>

namespace nonlocal::solver_2d::thermal {

template<std::floating_point T, std::integral I>
problem_settings init_problem_settings(const mesh::mesh_container_2d<T, I>& mesh,
                                       const parameters_2d<T>& parameters,
                                       const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                       const bool is_stationary_problem) {
    static constexpr auto is_flux = [](const auto& condition) noexcept {
        return  bool(dynamic_cast<const flux_2d<T>*>(condition.get())) &&
               !bool(dynamic_cast<const combined_flux_2d<T>*>(condition.get()));
    };
    static constexpr auto is_radiation_boundary = [](const auto& condition) noexcept {
        return bool(dynamic_cast<const radiation_2d<T>*>(condition.get()));
    };
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) {
        static constexpr auto is_not_constant = [](const coefficient_t<T, 2>& coefficient) { return !is_constant<T, 2>(coefficient); };
        return std::visit(metamath::types::visitor{is_not_constant,
            [](const auto& conductivity) { return std::any_of(conductivity.begin(), conductivity.end(), is_not_constant); }
        }, parameter.physical.conductivity);
    };
    static constexpr auto is_solution_dependent = [](const auto& parameter) {
        static constexpr auto is_solution_dependent = [](const coefficient_t<T, 2>& coefficient) { 
            return std::holds_alternative<solution_dependency<T, 2>>(coefficient); 
        };
        return std::visit(metamath::types::visitor{is_solution_dependent,
            [](const auto& conductivity) { return std::any_of(conductivity.begin(), conductivity.end(), is_solution_dependent); },
        }, parameter.physical.conductivity);
    };
    const auto parameters_view = parameters | std::views::values;
    const auto boundaries_view = boundaries_conditions | std::views::values;
    return {
        .theories = theories_types(parameters),
        .is_neumann = is_stationary_problem ? std::all_of(boundaries_view.begin(), boundaries_view.end(), is_flux) : false,
        .is_nonlinear_boundary = std::any_of(boundaries_view.begin(), boundaries_view.end(), is_radiation_boundary),
        .is_nonconstant_parameters = std::any_of(parameters_view.begin(), parameters_view.end(), is_nonconstant_parameters),
        .is_solution_dependent = std::any_of(parameters_view.begin(), parameters_view.end(), is_solution_dependent),
        .is_inner_nodes = utils::inner_nodes(mesh, boundaries_conditions)
    };
}

}
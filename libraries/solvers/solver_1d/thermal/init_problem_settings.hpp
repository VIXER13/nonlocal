#pragma once

#include "thermal_boundary_conditions_1d.hpp"
#include "thermal_parameters_1d.hpp"

#include <solvers/solver_1d/base/problem_settings.hpp>

namespace nonlocal::solver_1d::thermal {

template<class T>
problem_settings init_problem_settings(const parameters_1d<T>& parameters,
                                       const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                       const bool is_stationary_problem) {
    static constexpr auto is_flux = [](const auto& condition) noexcept {
        return  bool(dynamic_cast<const flux_1d<T>*>(condition.get())) &&
               !bool(dynamic_cast<const combined_flux_1d<T>*>(condition.get()));
    };
    static constexpr auto is_radiation_boundary = [](const auto& condition) noexcept {
        return bool(dynamic_cast<const radiation_1d<T>*>(condition.get()));
    };
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) noexcept {
        return std::holds_alternative<spatial_dependency<T, 1>>(parameter.physical.conductivity) ||
               std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.conductivity);
    };
    static constexpr auto is_solution_dependent = [](const auto& parameter) noexcept {
        return std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.conductivity);
    };
    return {
        .theories = theories_types(parameters),
        .is_neumann = is_stationary_problem ? std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_flux) : false,
        .is_nonconstant_parameters = std::any_of(parameters.begin(), parameters.end(), is_nonconstant_parameters),
        .is_radiation_boundary = std::any_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_radiation_boundary),
        .is_solution_dependent = std::any_of(parameters.begin(), parameters.end(), is_solution_dependent),
        .is_first_kind = { 
            bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.front().get())),
            bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.back ().get()))
        }
    };
}

}
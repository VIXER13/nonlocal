#pragma once

#include "boundary_data_1d.hpp"

#include <constants/nonlocal_constants.hpp>
#include <logger/logger.hpp>

#include <algorithm>
#include <array>
#include <vector>

namespace nonlocal::solver_1d {

struct problem_settings final {
    std::vector<theory_t> theories;
    bool is_neumann = false;
    bool is_nonlinear_boundary = false;
    bool is_nonconstant_parameters = false;
    bool is_solution_dependent = false;
    std::array<bool, 2> is_first_kind = { false, false };

    constexpr bool is_nonlinear() const noexcept {
        return is_nonlinear_boundary || is_solution_dependent;
    }

    bool is_symmetric() const {
        static constexpr auto is_nonlocal = [](const theory_t theory) noexcept { return theory == theory_t::NONLOCAL; };
        const bool is_any_nonlocal = std::any_of(theories.begin(), theories.end(), is_nonlocal);
        return !(is_any_nonlocal && is_nonconstant_parameters);
    }
};

template<class T>
problem_settings init_boundary_settings(problem_settings settings, const boundaries_data_1d<T>& boundaries,
                                         const bool is_stationary_problem) {
    static constexpr auto is_pure_load = [](const auto& condition) constexpr noexcept {
        return condition.is_pure_load();
    };
    static constexpr auto is_nonlinear = [](const auto& condition) constexpr noexcept {
        return condition.nonlinear.has_value();
    };
    settings.is_neumann = is_stationary_problem && std::all_of(boundaries.begin(), boundaries.end(), is_pure_load);
    settings.is_nonlinear_boundary = std::any_of(boundaries.begin(), boundaries.end(), is_nonlinear);
    settings.is_first_kind = { boundaries.front().prescribed_value.has_value(), boundaries.back().prescribed_value.has_value() };
    return settings;
}

inline void log_problem_settings(const problem_settings& settings) {
    logger::info() << (settings.is_symmetric() ? "Symmetric problem" : "Asymmetrical problem") << std::endl;
    if (settings.is_neumann)
        logger::info() << "Neuman problem" << std::endl;
    if (settings.is_nonlinear())
        logger::info() << "Nonlinear problem" << std::endl;
    if (settings.is_solution_dependent)
        throw std::domain_error{ "Parametrically nonlinear problems are not supported at the moment." };
}

} // namespace nonlocal::solver_1d

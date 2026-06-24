#pragma once

#include <constants/nonlocal_constants.hpp>
#include <logger/logger.hpp>

#include <algorithm>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

namespace nonlocal::solver_2d {

struct problem_settings final {
    std::unordered_map<std::string, theory_t> theories;
    bool is_neumann = false;
    bool is_nonlinear_boundary = false;
    bool is_nonconstant_parameters = false;
    bool is_solution_dependent = false;
    std::vector<bool> is_inner_nodes;

    constexpr bool is_nonlinear() const noexcept {
        return is_nonlinear_boundary || is_solution_dependent;
    }

    bool is_nonlocal() const {
        static constexpr auto is_nonlocal_check = [](const theory_t theory) noexcept { return theory == theory_t::NONLOCAL; };
        const auto theories_view = theories | std::views::values;
        return std::any_of(theories_view.begin(), theories_view.end(), is_nonlocal_check);
    }

    bool is_symmetric() const {
        return !(is_nonlocal() && is_nonconstant_parameters);
    }
};

inline void log_problem_settings(const problem_settings& settings) {
    logger::info() << (settings.is_symmetric() ? "Symmetric problem" : "Asymmetrical problem") << std::endl;
    if (settings.is_neumann)
        logger::info() << "Neuman problem" << std::endl;
    if (settings.is_nonlinear())
        logger::info() << "Nonlinear problem" << std::endl;
    if (settings.is_solution_dependent)
        throw std::domain_error{"Parametrically nonlinear problems are not supported at the moment."};
}

}
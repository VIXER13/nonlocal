#pragma once

#include <constants/nonlocal_constants.hpp>
#include <logger/logger.hpp>

#include <array>
#include <algorithm>
#include <vector>

namespace nonlocal {

struct problem_settings final {
    std::vector<theory_t> theories;
    bool is_neumann = false;
    bool is_nonconstant_parameters = false;
    bool is_radiation_boundary = false;
    bool is_solution_dependent = false;
    std::array<bool, 2> is_first_kind = {false, false};

    constexpr bool is_nonlinear() const noexcept {
        return is_radiation_boundary || is_solution_dependent;
    }

    bool is_symmetric() const {
        static constexpr auto is_nonlocal = [](const theory_t theory) noexcept {
            return theory == theory_t::NONLOCAL;
        };
        const bool is_any_nonlocal = std::any_of(theories.begin(), theories.end(), is_nonlocal);
        return !(is_any_nonlocal && is_nonconstant_parameters);
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
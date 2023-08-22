#ifndef NONLOCFEM_TWO_DIMENSIONAM_PROBLEMS_HPP
#define NONLOCFEM_TWO_DIMENSIONAM_PROBLEMS_HPP

#include "init_utils.hpp"

namespace nonlocal {

template<class T, class I>
void two_dimensional_problems(const nlohmann::json& config, const nonlocal::config::save_data& save, const config::problem_t problem) {
    static const std::set<config::problem_t> available_problems = {
        config::problem_t::THERMAL_STATIONARY,
        config::problem_t::THERMAL_NONSTATIONARY,
        config::problem_t::MECHANICAL_EQUILIBRIUM
    };
    if (!available_problems.contains(problem)) {
        throw std::domain_error{
            "In the two-dimensional case, the following problems are available: " +
            init_available_problems_list(available_problems)
        };
    }

    if (problem == config::problem_t::THERMAL_STATIONARY)
        std::cout << "thermal_stationary_2d" << std::endl;
    else if (problem == config::problem_t::THERMAL_NONSTATIONARY)
        std::cout << "thermal_nonstationary_2d" << std::endl;
    else if (problem == config::problem_t::MECHANICAL_EQUILIBRIUM)
        std::cout << "equilibrium" << std::endl;
}

}

#endif
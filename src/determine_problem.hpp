#ifndef NONLOCFEM_DETERMINE_PROBLEM_HPP
#define NONLOCFEM_DETERMINE_PROBLEM_HPP

#include "thermal_problems_1d.hpp"
#include "thermal_problems_2d.hpp"

#include <set>

namespace nonlocal {

class _determine_problem final {
    constexpr explicit _determine_problem() noexcept = default;

    static void init_save_data(const config::save_data& save, const nlohmann::json& config);

    static std::string init_available_problems_list(const std::set<config::problem_t>& available_problems);
    static std::vector<std::string> get_required_fields(const uint64_t dimension, const bool is_time_dependent);

    static bool is_thermal_problem(const config::problem_t problem);
    static bool is_mechanical_problem(const config::problem_t problem);
    static bool is_time_dependent_problem(const config::problem_t problem);

public:
    template<std::floating_point T, std::signed_integral I>
    friend void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem);

    template<std::floating_point T, std::signed_integral I>
    friend void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem);

    template<std::floating_point T, std::signed_integral I>
    friend void determine_problem(const nlohmann::json& config);
};

template<std::floating_point T, std::signed_integral I>
void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem) {
    static const std::set<config::problem_t> available_problems = {
        config::problem_t::THERMAL_STATIONARY,
        config::problem_t::THERMAL_NONSTATIONARY
    };
    if (!available_problems.contains(problem)) {
        throw std::domain_error{
            "In the one-dimensional case, the following problems are available: " +
            _determine_problem::init_available_problems_list(available_problems)
        };
    }

    config::check_required_fields(config, _determine_problem::get_required_fields(1, _determine_problem::is_time_dependent_problem(problem)));
    config::check_optional_fields(config, {"mesh", "auxiliary"});
    if (_determine_problem::is_thermal_problem(problem))
        thermal::solve_thermal_1d_problem<T, I>(config, save, problem);
}

template<std::floating_point T, std::signed_integral I>
void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem) {
    static const std::set<config::problem_t> available_problems = {
        config::problem_t::THERMAL_STATIONARY,
        config::problem_t::THERMAL_NONSTATIONARY,
        config::problem_t::MECHANICAL_EQUILIBRIUM
    };
    if (!available_problems.contains(problem)) {
        throw std::domain_error{
            "In the two-dimensional case, the following problems are available: " +
            _determine_problem::init_available_problems_list(available_problems)
        };
    }

    config::check_required_fields(config, _determine_problem::get_required_fields(2, _determine_problem::is_time_dependent_problem(problem)));
    config::check_optional_fields(config, {"auxiliary"});
    if (_determine_problem::is_thermal_problem(problem))
        thermal::solve_thermal_2d_problem<T, I>(config, save, problem);
}

template<std::floating_point T, std::signed_integral I>
void determine_problem(const nlohmann::json& config) {
    const bool contains_save = config.contains("save");
    const auto save = contains_save ? config::save_data{config["save"], "save"} : nonlocal::config::save_data{};
    if (contains_save)
        _determine_problem::init_save_data(save, config);
    else
        std::cerr << "WARNING: There is no \"save\" field in the config. Required data may not be saved." << std::endl;

    config::check_required_fields(config, {"task"});
    if (const config::task_data task{config["task"], "task"}; task.dimension == 1)
        problems_1d<T, I>(config, save, task.problem);
    else if (task.dimension == 2)
        problems_2d<T, I>(config, save, task.problem);
    else throw std::domain_error{"The dimension of the problem cannot be equal to " + std::to_string(task.dimension)};
}

}

#endif
#ifndef NONLOCFEM_DETERMINE_PROBLEM_HPP
#define NONLOCFEM_DETERMINE_PROBLEM_HPP

#include "thermal_problems_1d.hpp"
#include "thermal_problems_2d.hpp"
#include "mechanical_problems_2d.hpp"

#include <set>

namespace nonlocal {

class _determine_problem final {
    constexpr explicit _determine_problem() noexcept = default;

    static void init_save_data(const config::save_data& save, const nlohmann::json& config);

public:
    template<std::floating_point T, std::signed_integral I>
    friend void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T, std::signed_integral I>
    friend void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T, std::signed_integral I>
    friend void determine_problem(const nlohmann::json& config);
};

template<std::floating_point T, std::signed_integral I>
void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    config::check_required_fields(config, {"boundaries", "materials"});
    config::check_optional_fields(config, {"mesh", "auxiliary"});
    if (task.problem == nonlocal::config::problem_t::THERMAL)
        thermal::solve_thermal_1d_problem<T, I>(config, save, task.time_dependency);
    else throw std::domain_error{"Unknown task. In the one-dimensional case, the following problems are available: \"thermal\""};
}

template<std::floating_point T, std::signed_integral I>
void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    config::check_required_fields(config, {"boundaries", "materials", "mesh"});
    config::check_optional_fields(config, {"auxiliary"});
    auto mesh = std::make_shared<mesh::mesh_2d<T, I>>(config::mesh_data<2>{config["mesh"], "mesh"}.path);
    if (task.problem == nonlocal::config::problem_t::THERMAL)
        thermal::solve_thermal_2d_problem(mesh, config, save, task.time_dependency);
    else if (task.problem == nonlocal::config::problem_t::MECHANICAL)
        mechanical::solve_mechanical_2d_problem(mesh, config, save, task.time_dependency);
    else throw std::domain_error{"Unknown task. In the two-dimensional case, the following problems are available: \"thermal\", \"mechanical\""};
}

template<std::floating_point T, std::signed_integral I>
void determine_problem(const nlohmann::json& config) {
    const bool contains_save = config.contains("save");
    const auto save = contains_save ? config::save_data{config["save"], "save"} : nonlocal::config::save_data{};
    if (contains_save)
        _determine_problem::init_save_data(save, config);
    else
        logger::get().log(logger::log_level::WARNING) << "There is no \"save\" field in the config. Required data may not be saved." << std::endl;

    config::check_required_fields(config, {"task"});
    if (const config::task_data task{config["task"], "task"}; task.dimension == 1)
        problems_1d<T, I>(config, save, task);
    else if (task.dimension == 2)
        problems_2d<T, I>(config, save, task);
    else throw std::domain_error{"The dimension of the problem cannot be equal to " + std::to_string(task.dimension)};
}

}

#endif
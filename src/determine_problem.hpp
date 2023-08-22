#ifndef NONLOCFEM_DETERMINE_PROBLEM_HPP
#define NONLOCFEM_DETERMINE_PROBLEM_HPP

#include "one_dimensional_problems.hpp"
#include "two_dimensional_problems.hpp"

namespace nonlocal {

template<class T, class I>
void determine_problem(const nlohmann::json& config) {
    const bool contains_save = config.contains("save");
    const auto save = contains_save ? nonlocal::config::save_data{config["save"], "save"} : nonlocal::config::save_data{};
    if (contains_save) init_save_data(save, config);
    else std::cerr << "WARNING: There is no \"save\" field in the config. Required data may not be saved." << std::endl;

    nonlocal::config::check_required_fields(config, {"task"});
    if (const nonlocal::config::task_data task{config["task"], "task"}; task.dimension == 1)
        one_dimensional_problems<T, I>(config, save, task.problem);
    else if (task.dimension == 2)
        two_dimensional_problems<T, I>(config, save, task.problem);
    else throw std::domain_error{"The dimension of the problem cannot be equal to " + std::to_string(task.dimension)};
}

}

#endif
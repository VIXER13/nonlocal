#include "task_data.hpp"

#include "config_utils.hpp"

namespace nonlocal::config {

task_data::task_data(const nlohmann::json& config, const std::string& path) {
    nonlocal::config::check_required_fields(config, {"dimension", "problem", "time_dependency"}, append_access_sign(path));
    dimension = config["dimension"].get<size_t>();
    problem = config["problem"].get<problem_t>();
    time_dependency = config["time_dependency"].get<bool>();
}

task_data::operator nlohmann::json() const {
    return {
        {"dimension", dimension},
        {"problem", problem},
        {"time_dependency", time_dependency}
    };
}

}
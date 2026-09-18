#include "task_data.hpp"

#include "config_utils.hpp"

namespace nonlocal::config {

task_data::task_data(const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, {"dimension", "problem", "analysis_type"}, append_access_sign(path));
    dimension = config["dimension"].get<size_t>();
    problem = config["problem"].get<problem_t>();
    analysis_type = config["analysis_type"].get<analysis_type_t>();
}

task_data::operator nlohmann::json() const {
    return {
        {"dimension", dimension},
        {"problem", problem},
        {"analysis_type", analysis_type}
    };
}

}
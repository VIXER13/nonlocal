#include "determine_problem.hpp"

namespace nonlocal {

void _determine_problem::init_save_data(const nonlocal::config::save_data& save, const nlohmann::json& config) {
    if (!std::filesystem::exists(save.folder()))
        std::filesystem::create_directories(save.folder());
    if (save.contains("config"))
        nonlocal::config::dump_json(config, save.path("config", "json"));
}

std::string _determine_problem::init_available_problems_list(const std::set<config::problem_t>& available_problems) {
    std::string result = "[ ";
    for(auto it = available_problems.cbegin(); it != available_problems.cend(); ++it) {
        if (it != available_problems.cbegin())
            result += ", ";
        result += nlohmann::json(*it).get<std::string>();
    }
    return result + " ]";
}

std::vector<std::string> _determine_problem::get_required_fields(const bool is_time_dependent) {
    if (is_time_dependent)
        return {"boundaries", "materials", "time"};
    else
        return {"boundaries", "materials"};
}

bool _determine_problem::is_thermal_problem(const config::problem_t problem) {
    return problem == config::problem_t::THERMAL_STATIONARY ||
           problem == config::problem_t::THERMAL_NONSTATIONARY;
}

bool _determine_problem::is_mechanical_problem(const config::problem_t problem) {
    return problem == config::problem_t::MECHANICAL_EQUILIBRIUM;
}

bool _determine_problem::is_time_dependent_problem(const config::problem_t problem) {
    return problem == config::problem_t::THERMAL_NONSTATIONARY;
}

}
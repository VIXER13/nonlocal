#include "init_utils.hpp"

namespace nonlocal {

void init_save_data(const nonlocal::config::save_data& save, const nlohmann::json& config) {
    if (!std::filesystem::exists(save.folder()))
        std::filesystem::create_directories(save.folder());
    if (save.contains("config"))
        nonlocal::config::dump_json(config, save.path("config", "json"));
}

std::string init_available_problems_list(const std::set<config::problem_t>& available_problems) {
    std::string result = "[ ";
    for(auto it = available_problems.cbegin(); it != available_problems.cend(); ++it) {
        if (it != available_problems.cbegin())
            result += ", ";
        result += nlohmann::json(*it).get<std::string>();
    }
    return result + " ]";
}

}
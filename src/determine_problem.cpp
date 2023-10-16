#include "determine_problem.hpp"

namespace nonlocal {

void _determine_problem::init_save_data(const config::save_data& save, const nlohmann::json& config) {
    if (!std::filesystem::exists(save.folder()))
        std::filesystem::create_directories(save.folder());
    if (save.contains("config"))
        config::dump_json(config, save.path("config", "json"));
}

}
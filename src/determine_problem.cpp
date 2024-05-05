#include "determine_problem.hpp"

namespace nonlocal {

void _determine_problem::init_save_data(const config::save_data& save, const nlohmann::json& config) {
    if (!std::filesystem::exists(save.folder()))
        std::filesystem::create_directories(save.folder());
    if (save.contains("config"))
        config::dump_json(config, save.path("config", "json"));
}

std::vector<std::string> _determine_problem::get_required_fields(const config::task_data& task) {
    auto required_fields = task.problem == nonlocal::config::problem_t::THERMOMECHANICAL ? 
        std::vector<std::string>{"thermal_boundaries", "mechanical_boundaries", "thermal_materials", "mechanical_materials", "mesh"} :
        std::vector<std::string>{"boundaries", "materials", "mesh"};
    if (task.time_dependency)
        required_fields.push_back("time");
    return required_fields;
}

bool _determine_problem::is_thermal(const config::problem_t problem) {
    return problem == config::problem_t::THERMAL ||
           problem == config::problem_t::THERMOMECHANICAL;
}

bool _determine_problem::is_mechanical(const config::problem_t problem) {
    return problem == config::problem_t::MECHANICAL ||
           problem == config::problem_t::THERMOMECHANICAL;
}

}
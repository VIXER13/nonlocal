#include "thermal_config_data.hpp"

#include <iostream>

namespace {

enum class problem_t : uint8_t {
    UNKNOWN,
    THERMAL_STATIONARY,
    THERMAL_NONSTATIONARY,
    MECHANICAL_EQUILIBRIUM
};

NLOHMANN_JSON_SERIALIZE_ENUM(problem_t, {
    {problem_t::UNKNOWN, nullptr},
    {problem_t::THERMAL_STATIONARY, "thermal_stationary"},
    {problem_t::THERMAL_NONSTATIONARY, "thermal_nonstationary"},
    {problem_t::MECHANICAL_EQUILIBRIUM, "equilibrium"},
})

struct task_info final {
    problem_t problem = problem_t::UNKNOWN;
    uint8_t dimension = 0;

    explicit task_info(const nlohmann::json& config) {
        nonlocal::config::check_required_fields(config, {"task"});
        const nlohmann::json& task = config["task"];
        nonlocal::config::check_required_fields(task, {"problem", "dimension"}, "task");
        dimension = task["dimension"].get<uint>();
        problem = task["problem"].get<problem_t>();
    }
};

template<class T, class I>
void one_dimensional_problems(const nlohmann::json& config, const problem_t problem) {
    switch (problem)
    {
    case problem_t::THERMAL_STATIONARY:
        std::cout << "thermal_stationary_1d" << std::endl;
    break;
    
    case problem_t::THERMAL_NONSTATIONARY:
        std::cout << "thermal_nonstationary_1d" << std::endl;
    break;

    default:
        throw std::domain_error{
            "In the one-dimensional case, the following problems are available: "
            "[ " +
                nlohmann::json(problem_t::THERMAL_STATIONARY).get<std::string>() + ", " +
                nlohmann::json(problem_t::THERMAL_NONSTATIONARY).get<std::string>() +
            " ]"
        };
    }
}

template<class T, class I>
void two_dimensional_problems(const nlohmann::json& config, const problem_t problem) {
    switch (problem)
    {
    case problem_t::THERMAL_STATIONARY:
        std::cout << "thermal_stationary_2d" << std::endl;
    break;
    
    case problem_t::THERMAL_NONSTATIONARY:
        std::cout << "thermal_nonstationary_2d" << std::endl;
    break;

    case problem_t::MECHANICAL_EQUILIBRIUM:
        std::cout << "equilibrium_2d" << std::endl;
    break;

    default:
        throw std::domain_error{
            "In the one-dimensional case, the following problems are available: "
            "[ " +
                nlohmann::json(problem_t::THERMAL_STATIONARY).get<std::string>() + ", " +
                nlohmann::json(problem_t::THERMAL_NONSTATIONARY).get<std::string>() + ", " +
                nlohmann::json(problem_t::MECHANICAL_EQUILIBRIUM).get<std::string>() +
            " ]"
        };
    }
}

}

int main(const int argc, const char *const *const argv) {
    if (argc != 2) {
        std::cerr << "Input format: [program name] path/to/config.json" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        using T = double;
        using I = int64_t;
        const nlohmann::json config = nonlocal::config::parse_json(std::filesystem::path{argv[1]});
        switch (const task_info task{config}; task.dimension) {
        case 1:
            one_dimensional_problems<T, I>(config, task.problem);
        break;

        case 2:
            two_dimensional_problems<T, I>(config, task.problem);
        break;
        
        default:
            throw std::domain_error{"The dimension of the problem cannot be equal to " + std::to_string(task.dimension)};
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
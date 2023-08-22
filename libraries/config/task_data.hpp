#ifndef NONLOCAL_CONFIG_PROBLEM_DATA_HPP
#define NONLOCAL_CONFIG_PROBLEM_DATA_HPP

#include <nlohmann/json.hpp>

namespace nonlocal::config {

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

struct task_data final {
    uint64_t dimension = 0;
    problem_t problem = problem_t::UNKNOWN;

    explicit task_data(const nlohmann::json& config, const std::string& path);

    operator nlohmann::json() const;
};

}

#endif
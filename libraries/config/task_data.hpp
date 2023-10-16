#ifndef NONLOCAL_CONFIG_PROBLEM_DATA_HPP
#define NONLOCAL_CONFIG_PROBLEM_DATA_HPP

#include <nlohmann/json.hpp>

namespace nonlocal::config {

enum class problem_t : uint8_t {
    UNKNOWN,
    THERMAL,
    MECHANICAL,
    THERMOMECHANICAL
};

NLOHMANN_JSON_SERIALIZE_ENUM(problem_t, {
    {problem_t::UNKNOWN, nullptr},
    {problem_t::THERMAL, "thermal"},
    {problem_t::MECHANICAL, "mechanical"},
    {problem_t::THERMOMECHANICAL, "thermomechanical"},
})

struct task_data final {
    uint64_t dimension = 0;
    problem_t problem = problem_t::UNKNOWN;
    bool time_dependency = false;

    explicit task_data(const nlohmann::json& config, const std::string& path = {});

    operator nlohmann::json() const;
};

}

#endif
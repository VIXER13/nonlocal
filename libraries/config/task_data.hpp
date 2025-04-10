#pragma once

#include <nlohmann/json.hpp>

namespace nonlocal::config {

enum class problem_t : uint8_t {
    Unknown,
    Thermal,
    Mechanical,
    Thermomechanical
};

NLOHMANN_JSON_SERIALIZE_ENUM(problem_t, {
    {problem_t::Unknown, nullptr},
    {problem_t::Thermal, "thermal"},
    {problem_t::Mechanical, "mechanical"},
    {problem_t::Thermomechanical, "thermomechanical"},
})

struct task_data final {
    size_t dimension = 0;
    problem_t problem = problem_t::Unknown;
    bool time_dependency = false;

    explicit task_data(const nlohmann::json& config, const std::string& path = {});

    operator nlohmann::json() const;
};

}
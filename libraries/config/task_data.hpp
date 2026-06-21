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

enum class analysis_type_t : uint8_t {
    Stationary,
    TimeHarmonic,
    TimeDependent,
    Unknown
};

NLOHMANN_JSON_SERIALIZE_ENUM(analysis_type_t, {
    {analysis_type_t::Stationary, "stationary"},
    {analysis_type_t::TimeHarmonic, "time_harmonic"},
    {analysis_type_t::TimeDependent, "time_dependent"},
    {analysis_type_t::Unknown, nullptr}
})

struct task_data final {
    size_t dimension = 0;
    problem_t problem = problem_t::Unknown;
    analysis_type_t analysis_type = analysis_type_t::Unknown;

    explicit task_data(const nlohmann::json& config, const std::string& path = {});

    operator nlohmann::json() const;
};

}
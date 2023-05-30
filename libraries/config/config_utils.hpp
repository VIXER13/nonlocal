#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include "nonlocal_constants.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::config {

nlohmann::json parse_json(const std::filesystem::path& path);
void dump_json(const nlohmann::json& value, const std::filesystem::path& path);

// Throws an exception if at least one required field is missing.
void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required);

thermal::boundary_condition_t get_thermal_condition(const nlohmann::json& kind);
const std::string& get_thermal_condition(const thermal::boundary_condition_t kind);

size_t get_order(const nlohmann::json& order);
const std::string& get_order(const size_t order);
    
}

#endif
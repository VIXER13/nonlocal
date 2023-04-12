#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include "nonlocal_constants.hpp"

#include <json/value.h>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::config {

thermal::boundary_condition_t get_thermal_condition(const Json::Value& kind);
const std::string& get_thermal_condition(const thermal::boundary_condition_t kind);

size_t get_order(const Json::Value& order);
const std::string& get_order(const size_t order);

material_t get_material(const Json::Value& material);
const std::string& get_material(const material_t material);

Json::Value read_json(const std::filesystem::path& path);
Json::Value read_json(const char *const begin, const char *const end);
Json::Value read_json(const std::string& str);

void save_json(const std::filesystem::path& path, const Json::Value& value);

// Throws an exception if at least one required field is missing.
void check_required_fields(const Json::Value& value, const std::vector<std::string>& required);
    
}

#endif
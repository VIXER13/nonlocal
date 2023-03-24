#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include <json/json.h>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::config {

Json::Value read_json(const std::filesystem::path& config);
Json::Value read_json(const std::string& config);

// Throws an exception if at least one required field is missing.
void check_required_fields(const Json::Value& value, const std::vector<std::string>& required);
    
}

#endif
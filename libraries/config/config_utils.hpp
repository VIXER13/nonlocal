#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include <json/json.h>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::config {

Json::Value read_json(const std::filesystem::path& path);
Json::Value read_json(const char *const begin, const char *const end);
Json::Value read_json(const std::string& str);

void save_json(const std::filesystem::path& path, const Json::Value& value);

// Throws an exception if at least one required field is missing.
void check_required_fields(const Json::Value& value, const std::vector<std::string>& required);
    
}

#endif
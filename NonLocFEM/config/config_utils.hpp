#pragma once

#include <nlohmann/json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::config {

nlohmann::json parse_json(const std::filesystem::path& path);

void dump_json(const nlohmann::json& value, const std::filesystem::path& path, const int indent = 4, const char indent_char = ' ');

std::string append_access_sign(std::string path, const std::optional<size_t> index = std::nullopt);

// Throws an exception if at least one required field is missing.
// Specify a path for nested fields to get more descriptive error messages.
void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required, const std::string& path = {});

// Log info about missed optional fields
void check_optional_fields(const nlohmann::json& value, const std::vector<std::string>& optional, const std::string& path = {});
    
}
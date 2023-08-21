#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include "nonlocal_constants.hpp"

#include <nlohmann/json.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace nonlocal::thermal {

NLOHMANN_JSON_SERIALIZE_ENUM(boundary_condition_t, {
    {boundary_condition_t::TEMPERATURE, "temperature"},
    {boundary_condition_t::FLUX, "flux"},
    {boundary_condition_t::CONVECTION, "convection"},
    {boundary_condition_t::RADIATION, "radiation"},
    {boundary_condition_t::COMBINED, "combined"}
})

}

namespace nonlocal::mechanical {

NLOHMANN_JSON_SERIALIZE_ENUM(boundary_condition_t, {
    {boundary_condition_t::DISPLACEMENT, "displacement"},
    {boundary_condition_t::PRESSURE, "pressure"}
})

}

namespace nonlocal::config {

nlohmann::json parse_json(const std::filesystem::path& path);
void dump_json(const nlohmann::json& value, const std::filesystem::path& path, const int indent = 4, const char indent_char = ' ');

// Throws an exception if at least one required field is missing.
void check_required_fields(const nlohmann::json& value, const std::vector<std::string>& required, const std::string& path = "");

size_t get_order(const nlohmann::json& order);
const std::string& get_order(const size_t order);
    
}

#endif
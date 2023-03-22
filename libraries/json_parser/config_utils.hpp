#ifndef NONLOCAL_CONFIG_UTILS_HPP
#define NONLOCAL_CONFIG_UTILS_HPP

#include <json/json.h>

#include <filesystem>
#include <unordered_map>
#include <optional>
#include <string>
#include <vector>

namespace nonlocal::config {

class save_data final {
    std::filesystem::path _folder;
    std::unordered_map<std::string, std::string> _names;

public:
    explicit save_data(const Json::Value& save);

    bool contains(const std::string& key) const;
    const std::filesystem::path& folder() const noexcept;
    std::filesystem::path path(const std::string& key, const std::string& extension = ".csv",
                               const std::optional<std::string>& default_name = std::nullopt) const;
};

Json::Value read_json_from_file(const std::filesystem::path& config);

// Throws an exception if at least one required field is missing.
void check_required_fields(const Json::Value& value, const std::vector<std::string>& required);
    
}

#endif
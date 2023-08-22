#ifndef NONLOCAL_CONFIG_SAVE_DATA_HPP
#define NONLOCAL_CONFIG_SAVE_DATA_HPP

#include <nlohmann/json.hpp>

#include <filesystem>
#include <unordered_map>
#include <string>
#include <optional>

namespace nonlocal::config {

class save_data final {
    std::filesystem::path _folder = std::filesystem::current_path();
    std::unordered_map<std::string, std::string> _names;
    std::optional<std::streamsize> _precision;

public:
    explicit save_data() = default;
    explicit save_data(const nlohmann::json& config, const std::string& path = "");

    std::optional<std::streamsize> precision() const noexcept;
    const std::filesystem::path& folder() const noexcept;

    bool contains(const std::string& key) const;
    std::string get_name(const std::string& key, const std::optional<std::string>& default_name = std::nullopt) const;

    std::filesystem::path make_path(const std::string& name, const std::string& extension) const;
    std::filesystem::path path(const std::string& key, const std::string& extension,
                               const std::optional<std::string>& default_name = std::nullopt) const;

    operator nlohmann::json() const;
};

}

#endif
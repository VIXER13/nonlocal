#pragma once

#include "config_utils.hpp"

namespace nonlocal::config {

template<size_t Dimension>
struct mesh_data final {
    std::filesystem::path path; // required

    explicit mesh_data() = default;
    explicit mesh_data(const nlohmann::json& config, const std::string& config_path = {}) {
        check_required_fields(config, { "path" }, append_access_sign(config_path));
        path = config["path"].get<std::string>();
    }

    operator nlohmann::json() const {
        return { {"path", path.string()} };
    }
};

}
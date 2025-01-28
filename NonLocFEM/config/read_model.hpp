#pragma once

#include "config_utils.hpp"
#include "../equation_parameters.hpp"
#include "influence_functions_1d.hpp"

namespace nonlocal::config {

template<std::floating_point T>
model_parameters<1u, T> read_model_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, { "local_weight", "nonlocal_radius" }, path_with_access);
    return {
        .influence = influence::polynomial_1d<T, 1u, 1u>{config["nonlocal_radius"].get<T>()},
        .local_weight = config["local_weight"].get<T>()
    };
}

}
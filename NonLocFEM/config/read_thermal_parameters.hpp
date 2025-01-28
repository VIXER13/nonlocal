#pragma once

#include "read_model.hpp"

#include "thermal/thermal_parameters_1d.hpp"

namespace nonlocal::config {

template<std::floating_point T>
thermal::parameter_1d_sptr<T> read_thermal_coefficient_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"conductivity"}, path_with_access);
    check_optional_fields(config, {"capacity", "density", "relaxation_time"}, path_with_access);
    return std::make_shared<thermal::parameter_1d<T, coefficients_t::CONSTANTS>>(
        config["conductivity"].get<T>(),
        config.value("capacity", T{1}),
        config.value("density", T{1}),
        config.value("relaxation_time", T{0})
    );
}

template<std::floating_point T>
thermal::parameters_1d<T> read_thermal_parameters_1d(const nlohmann::json& config, const std::string& path) {
    if (config.empty() || !config.is_array())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an nonempty array."};

    thermal::parameters_1d<T> parameters(config.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.capacity()}) {
        const nlohmann::json& config_material = config[i];
        const std::string path_with_access = append_access_sign(append_access_sign(path, i));
        check_required_fields(config_material, {"physical"}, path_with_access);
        const std::string model_field = config_material.contains("thermal_model") ? "thermal_model" : 
                                        config_material.contains("model")         ?         "model" : "";
        if (model_field.empty())
            logger::debug() << "Optional field \"" + path_with_access + "model\" is not specified." << std::endl;
        parameters[i] = {
            .model = model_field.empty() ? model_parameters<1u, T>{} : read_model_1d<T>(config_material[model_field], path_with_access + model_field),
            .physical = read_thermal_coefficient_1d<T>(config_material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

}
#pragma once

#include "read_model.hpp"

#include "thermal/thermal_parameters_1d.hpp"
#include "thermal/thermal_parameters_2d.hpp"

namespace nonlocal::config {

template<std::floating_point T>
thermal::parameter_1d_sptr<T> read_thermal_coefficient_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"conductivity"}, path_with_access);
    check_optional_fields(config, {"capacity", "density", "relaxation_time"}, path_with_access);
    return std::make_shared<thermal::parameter_1d<T>>(
        config["conductivity"].get<T>(),
        config.value("capacity", T{1}),
        config.value("density", T{1}),
        config.value("relaxation_time", T{0})
    );
}

template<std::floating_point T>
thermal::parameters_1d<T> read_thermal_parameters_1d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_array())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an nonempty array."};

    thermal::parameters_1d<T> parameters(config.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.capacity()}) {
        const nlohmann::json& config_material = config[i];
        const std::string path_with_access = append_access_sign(append_access_sign(path, i));
        check_required_fields(config_material, {"physical"}, path_with_access);
        const std::string model_field = get_model_field(config_material, path_with_access, "thermal");
        parameters[i] = {
            .model = model_field.empty() ? model_parameters<1u, T>{} : read_model_1d<T>(config_material[model_field], path_with_access + model_field),
            .physical = read_thermal_coefficient_1d<T>(config_material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

template<std::floating_point T>
metamath::types::square_matrix<T, 2> read_conductivity_2d(const nlohmann::json& config, const std::string& path) {
    if (config.is_number()) {
        const T value = config.get<T>();
        return {
            value, T{0},
            T{0}, value
        };
    }
    if (config.is_array() && config.size() == 2)
        return {
            config[0].get<T>(), T{0},
            T{0}, config[1].get<T>()
        };
    if (config.is_array() && config.size() == 4)
        return {
            config[0].get<T>(), config[1].get<T>(),
            config[2].get<T>(), config[3].get<T>()
        };
    throw std::domain_error{"The thermal conductivity should be either a number in the isotropic case, "
                            "or an array of size 2 in the orthotropic case and size 4 in the anisotropic case."};
}

template<std::floating_point T>
thermal::parameter_2d_sptr<T> read_thermal_coefficient_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"conductivity"}, path_with_access);
    check_optional_fields(config, {"capacity", "density", "relaxation_time"}, path_with_access);
    return std::make_shared<thermal::parameter_2d<T>>(
        read_conductivity_2d<T>(config["conductivity"], path_with_access + "conductivity"),
        config.value("capacity", T{1}),
        config.value("density", T{1})
    );
}

template<std::floating_point T>
thermal::parameters_2d<T> read_thermal_parameters_2d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_object())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an object."};
    const std::string path_with_access = append_access_sign(path);
    thermal::parameters_2d<T> parameters;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_access_to_material = append_access_sign(path_with_access + name);
        const std::string model_field = get_model_field(material, path_with_access, "thermal");
        auto& parameter = parameters[name] = {
            .model = model_field.empty() ? model_parameters<2u, T>{} : read_model_2d<T>(material[model_field], path_with_access + model_field),
            .physical = read_thermal_coefficient_2d<T>(material["physical"], path_with_access + "physical")
        };
        parameter.physical->material = material_t::ISOTROPIC;
    }
    return parameters;
}

}
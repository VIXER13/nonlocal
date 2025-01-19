#pragma once

#include "config_utils.hpp"

#include "thermal/thermal_boundary_conditions_1d.hpp"
#include "thermal/thermal_boundary_conditions_2d.hpp"

namespace nonlocal::config {

enum class thermal_boundary_condition_t : uint8_t {
    Undefined,
    Temperature,
    Flux,
    Convection,
    Radiation,
    Combined
};

NLOHMANN_JSON_SERIALIZE_ENUM(thermal_boundary_condition_t, {
    {thermal_boundary_condition_t::Undefined, nullptr},
    {thermal_boundary_condition_t::Temperature, "temperature"},
    {thermal_boundary_condition_t::Flux, "flux"},
    {thermal_boundary_condition_t::Convection, "convection"},
    {thermal_boundary_condition_t::Radiation, "radiation"},
    {thermal_boundary_condition_t::Combined, "combined"}
})

template<std::floating_point T>
std::unique_ptr<thermal::thermal_boundary_condition_1d<T>> read_thermal_boundary_condition_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"kind"}, path_with_access);
    switch (config["kind"].get<thermal_boundary_condition_t>()) {
    case thermal_boundary_condition_t::Temperature:
        check_required_fields(config, { "temperature" }, path_with_access);
        return std::make_unique<thermal::temperature_1d<T>>(config["temperature"].get<T>());

    case thermal_boundary_condition_t::Flux:
        check_required_fields(config, { "flux" }, path_with_access);
        return std::make_unique<thermal::flux_1d<T>>(config["flux"].get<T>());

    case thermal_boundary_condition_t::Convection:
        check_required_fields(config, { "temperature", "heat_transfer" }, path_with_access);
        return std::make_unique<thermal::convection_1d<T>>(config["temperature"].get<T>(), config["heat_transfer"].get<T>());

    case thermal_boundary_condition_t::Radiation:
        check_required_fields(config, { "emissivity" }, path_with_access);
        return std::make_unique<thermal::radiation_1d<T>>(config["emissivity"].get<T>());

    case thermal_boundary_condition_t::Combined:
        if (!config.contains("heat_transfer"))
            check_optional_fields(config, {"temperature", "flux", "heat_transfer", "emissivity"}, path_with_access);
        else {
            check_required_fields(config, {"temperature"}, path_with_access);
            check_optional_fields(config, {"flux", "emissivity"}, path_with_access);
        }
        return std::make_unique<thermal::combined_flux_1d<T>>(
            config.value("flux", T{0}),
            config.value("heat_transfer", T{0}), config.value("temperature", T{0}),
            config.value("emissivity", T{0}));
    }
    throw std::domain_error{"Unknown boundary condition type: " + config["kind"].get<std::string>()};
}

template<std::floating_point T>
thermal::thermal_boundary_condition_2d<T> read_thermal_boundary_condition_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"kind"}, path_with_access);
    switch (config["kind"].get<thermal_boundary_condition_t>()) {
    case thermal_boundary_condition_t::Temperature:
        check_required_fields(config, { "temperature" }, path_with_access);
        return std::make_unique<thermal::temperature_2d<T>>(config["temperature"].get<T>());

    case thermal_boundary_condition_t::Flux:
        check_required_fields(config, { "flux" }, path_with_access);
        return std::make_unique<thermal::flux_2d<T>>(config["flux"].get<T>());

    case thermal_boundary_condition_t::Convection:
        check_required_fields(config, { "temperature", "heat_transfer" }, path_with_access);
        return std::make_unique<thermal::convection_2d<T>>(config["temperature"].get<T>(), config["heat_transfer"].get<T>());

    case thermal_boundary_condition_t::Radiation:
        check_required_fields(config, { "emissivity" }, path_with_access);
        return std::make_unique<thermal::radiation_2d<T>>(config["emissivity"].get<T>());

    case thermal_boundary_condition_t::Combined:
        if (!config.contains("heat_transfer"))
            check_optional_fields(config, {"temperature", "flux", "heat_transfer", "emissivity"}, path_with_access);
        else {
            check_required_fields(config, {"temperature"}, path_with_access);
            check_optional_fields(config, {"flux", "emissivity"}, path_with_access);
        }
        return std::make_unique<thermal::combined_flux_2d<T>>(
            config.value("flux", T{0}),
            config.value("heat_transfer", T{0}), config.value("temperature", T{0}),
            config.value("emissivity", T{0}));
    }
    throw std::domain_error{"Unknown boundary condition type: " + config["kind"].get<std::string>()};
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_1d<T> read_thermal_boundary_conditions_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"left", "right"}, path_with_access);
    return {
        read_thermal_boundary_condition_1d<T>(config["left"], path_with_access + "left"),
        read_thermal_boundary_condition_1d<T>(config["right"], path_with_access + "right")
    };
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_2d<T> read_thermal_boundary_conditions_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
    for(const auto& [name, condition] : config.items())
        boundaries_conditions[name] = read_thermal_boundary_condition_2d<T>(condition, path_with_access + name);
    return boundaries_conditions;
}

}
#pragma once

#include "config_utils.hpp"

#include <solvers/solver_1d/thermal/thermal_boundary_conditions_1d.hpp>
#include <solvers/solver_2d/thermal/thermal_boundary_conditions_2d.hpp>

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

class _read_thermal_boundary_conditions final {
    template<std::floating_point T>
    static void check_parameters(const T heat_transfer, const T emissivity, const std::string& path_with_access);

    template<std::floating_point T>
    static std::unique_ptr<thermal::thermal_boundary_condition_1d<T>> read_thermal_boundary_condition_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static thermal::thermal_boundary_condition_2d<T> read_thermal_boundary_condition_2d(const nlohmann::json& config, const std::string& path);

    explicit _read_thermal_boundary_conditions() noexcept = default;

public:
    template<std::floating_point T>
    friend thermal::thermal_boundaries_conditions_1d<T> read_thermal_boundaries_conditions_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend thermal::thermal_boundaries_conditions_2d<T> read_thermal_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
void _read_thermal_boundary_conditions::check_parameters(const T heat_transfer, const T emissivity, const std::string& path_with_access) {
    if (heat_transfer < T{0})
        throw std::domain_error{"\"" + path_with_access + "heat_transfer\" parameter shall be greather than 0."};
    if (emissivity < T{0} || emissivity > T{1})
        throw std::domain_error{"\"" + path_with_access + "emissivity\" parameter shall be in the interval [0, 1]."};
}

template<std::floating_point T>
std::unique_ptr<thermal::thermal_boundary_condition_1d<T>> 
_read_thermal_boundary_conditions::read_thermal_boundary_condition_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"kind"}, path_with_access);
    switch (config["kind"].get<thermal_boundary_condition_t>()) {
    case thermal_boundary_condition_t::Temperature:
        check_required_fields(config, { "temperature" }, path_with_access);
        return std::make_unique<thermal::temperature_1d<T>>(config["temperature"].get<T>());

    case thermal_boundary_condition_t::Flux:
        check_required_fields(config, { "flux" }, path_with_access);
        return std::make_unique<thermal::flux_1d<T>>(config["flux"].get<T>());

    case thermal_boundary_condition_t::Convection: {
        check_required_fields(config, { "temperature", "heat_transfer" }, path_with_access);
        const T heat_transfer = config["heat_transfer"].get<T>();
        static constexpr T emissivity = T{0};
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::convection_1d<T>>(config["temperature"].get<T>(), heat_transfer);
    }

    case thermal_boundary_condition_t::Radiation: {
        check_required_fields(config, { "emissivity" }, path_with_access);
        const T emissivity = config["emissivity"].get<T>();
        static constexpr T heat_transfer = T{0};
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::radiation_1d<T>>(emissivity);
    }

    case thermal_boundary_condition_t::Combined: {
        if (!config.contains("heat_transfer"))
            check_optional_fields(config, {"flux", "heat_transfer", "emissivity"}, path_with_access);
        else {
            check_required_fields(config, {"temperature"}, path_with_access);
            check_optional_fields(config, {"flux", "emissivity"}, path_with_access);
        }
        const T heat_transfer = config.value("heat_transfer", T{0});
        const T emissivity = config.value("emissivity", T{0});
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::combined_flux_1d<T>>(
            config.value("flux", T{0}),
            heat_transfer, config.value("temperature", T{0}),
            emissivity);
    }
    }
    throw std::domain_error{"Unknown boundary condition type: " + config["kind"].get<std::string>()};
}

template<std::floating_point T>
thermal::thermal_boundary_condition_2d<T> 
_read_thermal_boundary_conditions::read_thermal_boundary_condition_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"kind"}, path_with_access);
    switch (config["kind"].get<thermal_boundary_condition_t>()) {
    case thermal_boundary_condition_t::Temperature:
        check_required_fields(config, { "temperature" }, path_with_access);
        return std::make_unique<thermal::temperature_2d<T>>(config["temperature"].get<T>());

    case thermal_boundary_condition_t::Flux:
        check_required_fields(config, { "flux" }, path_with_access);
        return std::make_unique<thermal::flux_2d<T>>(config["flux"].get<T>());

    case thermal_boundary_condition_t::Convection: {
        check_required_fields(config, { "temperature", "heat_transfer" }, path_with_access);
        const T heat_transfer = config["heat_transfer"].get<T>();
        static constexpr T emissivity = T{0};
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::convection_2d<T>>(config["temperature"].get<T>(), heat_transfer);
    }

    case thermal_boundary_condition_t::Radiation: {
        check_required_fields(config, { "emissivity" }, path_with_access);
        const T emissivity = config["emissivity"].get<T>();
        static constexpr T heat_transfer = T{0};
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::radiation_2d<T>>(emissivity);
    }

    case thermal_boundary_condition_t::Combined: {
        if (!config.contains("heat_transfer"))
            check_optional_fields(config, {"flux", "heat_transfer", "emissivity"}, path_with_access);
        else {
            check_required_fields(config, {"temperature"}, path_with_access);
            check_optional_fields(config, {"flux", "emissivity"}, path_with_access);
        }
        const T heat_transfer = config.value("heat_transfer", T{0});
        const T emissivity = config.value("emissivity", T{0});
        check_parameters(heat_transfer, emissivity, path_with_access);
        return std::make_unique<thermal::combined_flux_2d<T>>(
            config.value("flux", T{0}),
            heat_transfer, config.value("temperature", T{0}),
            emissivity);
    }
    }
    throw std::domain_error{"Unknown boundary condition type: " + config["kind"].get<std::string>()};
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_1d<T> read_thermal_boundaries_conditions_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"left", "right"}, path_with_access);
    using _base = _read_thermal_boundary_conditions;
    return {
        _base::read_thermal_boundary_condition_1d<T>(config["left"], path_with_access + "left"),
        _base::read_thermal_boundary_condition_1d<T>(config["right"], path_with_access + "right")
    };
}

template<std::floating_point T>
thermal::thermal_boundaries_conditions_2d<T> read_thermal_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
    using _base = _read_thermal_boundary_conditions;
    for(const auto& [name, condition] : config.items())
        boundaries_conditions[name] = _base::read_thermal_boundary_condition_2d<T>(condition, path_with_access + name);
    return boundaries_conditions;
}

}
#ifndef NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP

#include "config_utils.hpp"

#include <exception>
#include <iostream>

namespace nonlocal::config {

enum class thermal_boundary_condition_t : uint8_t {
    UNDEFINED,
    TEMPERATURE,
    FLUX,
    CONVECTION,
    RADIATION,
    COMBINED
};

NLOHMANN_JSON_SERIALIZE_ENUM(thermal_boundary_condition_t, {
    {thermal_boundary_condition_t::UNDEFINED, nullptr},
    {thermal_boundary_condition_t::TEMPERATURE, "temperature"},
    {thermal_boundary_condition_t::FLUX, "flux"},
    {thermal_boundary_condition_t::CONVECTION, "convection"},
    {thermal_boundary_condition_t::RADIATION, "radiation"},
    {thermal_boundary_condition_t::COMBINED, "combined"}
})

template<std::floating_point T, size_t Dimension>
struct thermal_boundary_condition_data final {
    thermal_boundary_condition_t kind = thermal_boundary_condition_t::FLUX; // required
    T temperature = T{0};   // required if kind == TEMPERATURE or kind == CONVECTION
                            // used for TEMPERATURE condition, but if condition is CONVECTION used like ambient_temperature
    T flux = T{0};          // required if kind == FLUX
    T heat_transfer = T{0}; // required if kind == CONVECTION
    T emissivity = T{0};    // required if kind == RADIATION

    explicit constexpr thermal_boundary_condition_data() noexcept = default;
    explicit thermal_boundary_condition_data(const nlohmann::json& config, const std::string& path = {}) {
        const nlohmann::json& conf = config.contains("thermal") ? config["thermal"] : config;
        const std::string& config_path = config.contains("thermal") ? append_access_sign(path) + "thermal" : path;
        const std::string path_with_access = append_access_sign(config_path);
        check_required_fields(conf, { "kind" }, path_with_access);

        switch (kind = conf["kind"]) {
        case thermal_boundary_condition_t::TEMPERATURE: {
            check_required_fields(conf, { "temperature" }, path_with_access);
            temperature = conf["temperature"].get<T>();
        } break;

        case thermal_boundary_condition_t::FLUX: {
            check_required_fields(conf, { "flux" }, path_with_access);
            flux = conf["flux"].get<T>();
        } break;

        case thermal_boundary_condition_t::CONVECTION: {
            check_required_fields(conf, { "temperature", "heat_transfer" }, path_with_access);
            temperature = conf["temperature"].get<T>();
            heat_transfer = conf["heat_transfer"].get<T>();
        } break;

        case thermal_boundary_condition_t::RADIATION: {
            check_required_fields(conf, { "emissivity" }, path_with_access);
            emissivity = conf["emissivity"].get<T>();
        } break;

        case thermal_boundary_condition_t::COMBINED: {
            check_optional_fields(conf, {"temperature", "flux", "heat_transfer", "emissivity"}, path_with_access);
            temperature = conf.value("temperature", T{0});
            flux = conf.value("flux", T{0});
            heat_transfer = conf.value("heat_transfer", T{0});
            emissivity = conf.value("emissivity", T{0});
        } break;

        default:
            throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(kind))};
        }
    }

    operator nlohmann::json() const {
        return {
            {"kind", kind},
            {"temperature", temperature},
            {"flux", flux},
            {"heat_transfer", heat_transfer},
            {"emissivity", emissivity}
        };
    }
};

}

#endif
#ifndef NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP

#include "config_utils.hpp"

#include <exception>
#include <iostream>

namespace nonlocal::config {

enum class thermal_boundary_condition_t : uint8_t {
    UNKNOWN,
    TEMPERATURE,
    FLUX,
    CONVECTION,
    RADIATION,
    COMBINED
};

NLOHMANN_JSON_SERIALIZE_ENUM(thermal_boundary_condition_t, {
    {thermal_boundary_condition_t::UNKNOWN, nullptr},
    {thermal_boundary_condition_t::TEMPERATURE, "temperature"},
    {thermal_boundary_condition_t::FLUX, "flux"},
    {thermal_boundary_condition_t::CONVECTION, "convection"},
    {thermal_boundary_condition_t::RADIATION, "radiation"},
    {thermal_boundary_condition_t::COMBINED, "combined"}
})

template<std::floating_point T, size_t Dimension = 0> // Dimension == 0, because it does not depend on the dimension 
struct thermal_boundary_condition_data final {
    thermal_boundary_condition_t kind = thermal_boundary_condition_t::FLUX; // required
    T temperature = T{0};   // required if kind == TEMPERATURE or kind == CONVECTION
                            // used for TEMPERATURE condition, but if condition is CONVECTION used like ambient_temperature
    T flux = T{0};          // required if kind == FLUX
    T heat_transfer = T{0}; // required if kind == CONVECTION
    T emissivity = T{0};    // required if kind == RADIATION

    explicit constexpr thermal_boundary_condition_data() noexcept = default;
    explicit thermal_boundary_condition_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, { "kind" }, path_with_access);

        switch (kind = config["kind"]) {
        case thermal_boundary_condition_t::TEMPERATURE: {
            check_required_fields(config, { "temperature" }, path_with_access);
            temperature = config["temperature"].get<T>();
        } break;

        case thermal_boundary_condition_t::FLUX: {
            check_required_fields(config, { "flux" }, path_with_access);
            flux = config["flux"].get<T>();
        } break;

        case thermal_boundary_condition_t::CONVECTION: {
            check_required_fields(config, { "temperature", "heat_transfer" }, path_with_access);
            temperature = config["temperature"].get<T>();
            heat_transfer = config["heat_transfer"].get<T>();
        } break;

        case thermal_boundary_condition_t::RADIATION: {
            check_required_fields(config, { "emissivity" }, path_with_access);
            emissivity = config["emissivity"].get<T>();
        } break;

        case thermal_boundary_condition_t::COMBINED: {
            check_optional_fields(config, {"temperature", "flux", "heat_transfer", "emissivity"}, path_with_access);
            temperature = config.value("temperature", T{0});
            flux = config.value("flux", T{0});
            heat_transfer = config.value("heat_transfer", T{0});
            emissivity = config.value("emissivity", T{0});
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
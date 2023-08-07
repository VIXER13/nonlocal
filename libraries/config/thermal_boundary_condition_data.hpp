#ifndef NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP

#include "config_utils.hpp"

#include <exception>
#include <iostream>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension = 0>
struct thermal_boundary_condition_data final {
    thermal::boundary_condition_t kind = thermal::boundary_condition_t::FLUX; // required
    T temperature = T{0};   // required if kind == TEMPERATURE or kind == CONVECTION
                            // used for TEMPERATURE condition, but if condition is CONVECTION used like ambient_temperature
    T flux = T{0};          // required if kind == FLUX
    T heat_transfer = T{0}; // required if kind == CONVECTION
    T emissivity = T{0};    // required if kind == RADIATION

    explicit constexpr thermal_boundary_condition_data() noexcept = default;
    explicit thermal_boundary_condition_data(const nlohmann::json& condition) {
        check_required_fields(condition, { "kind" });

        switch (kind = condition["kind"]) {
        case thermal::boundary_condition_t::TEMPERATURE: {
            check_required_fields(condition, { "temperature" });
            temperature = condition["temperature"].get<T>();
        } break;

        case thermal::boundary_condition_t::FLUX: {
            check_required_fields(condition, { "flux" });
            flux = condition["flux"].get<T>();
        } break;

        case thermal::boundary_condition_t::CONVECTION: {
            check_required_fields(condition, { "temperature", "heat_transfer" });
            temperature = condition["temperature"].get<T>();
            heat_transfer = condition["heat_transfer"].get<T>();
        } break;

        case thermal::boundary_condition_t::RADIATION: {
            check_required_fields(condition, { "emissivity" });
            emissivity = condition["emissivity"].get<T>();
        } break;

        case thermal::boundary_condition_t::COMBINED: {
            temperature = condition.value("temperature", T{0});
            flux = condition.value("flux", T{0});
            heat_transfer = condition.value("heat_transfer", T{0});
            emissivity = condition.value("emissivity", T{0});
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
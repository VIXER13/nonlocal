#ifndef NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_BOUNDARY_CONDITION_DATA_HPP

#include <json/value.h>

#include <type_traits>
#include <exception>

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
    explicit thermal_boundary_condition_data(const Json::Value& condition) {
        check_required_fields(condition, { "kind" });
        switch (kind = get_thermal_condition(condition["kind"])) {
        case thermal::boundary_condition_t::TEMPERATURE: {
            check_required_fields(condition, { "temperature" });
            temperature = condition["temperature"].template as<T>();
        } break;

        case thermal::boundary_condition_t::FLUX: {
            check_required_fields(condition, { "flux" });
            flux = condition["flux"].template as<T>();
        } break;

        case thermal::boundary_condition_t::CONVECTION: {
            check_required_fields(condition, { "temperature", "heat_transfer" });
            temperature = condition["temperature"].template as<T>();
            heat_transfer = condition["heat_transfer"].template as<T>();
        } break;

        case thermal::boundary_condition_t::RADIATION: {
            check_required_fields(condition, { "emissivity" });
            emissivity = condition["emissivity"].template as<T>();
        } break;

        case thermal::boundary_condition_t::COMBINED: {
            temperature = condition.get("temperature", T{0}).template as<T>();
            flux = condition.get("flux", T{0}).template as<T>();
            heat_transfer = condition.get("heat_transfer", T{0}).template as<T>();
            emissivity = condition.get("emissivity", T{0}).template as<T>();
        } break;

        default:
            throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(kind))};
        }
    }

    Json::Value to_json() const {
        Json::Value result;
        result["kind"] = get_thermal_condition(kind);
        result["temperature"] = temperature;
        result["flux"] = flux;
        result["heat_transfer"] = heat_transfer;
        result["emissivity"] = emissivity;
        return result;
    }
};

}

#endif
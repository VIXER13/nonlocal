#ifndef NONLOCAL_THERMAL_CONFIG_DATA_HPP
#define NONLOCAL_THERMAL_CONFIG_DATA_HPP

#include "general_config_data.hpp"

#include "nonlocal_constants.hpp"

#include <exception>

namespace nonlocal::config {

thermal::boundary_condition_t get_thermal_condition(const Json::Value& kind);
size_t get_order(const Json::Value& order);

template<std::floating_point T, size_t Dimension>
struct thermal_material_data;

template<std::floating_point T>
struct thermal_material_data<T, 1> final {
    T conductivity = T{1}; // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const Json::Value& physical) {
        check_required_fields(physical, { "conductivity" });
        conductivity = physical["conductivity"].template as<T>();
        capacity = physical.get("capacity", T{1}).template as<T>();
        density = physical.get("density", T{1}).template as<T>();
    }
};

template<std::floating_point T>
struct thermal_equation_data final {
    T energy = 0;               // Used for Neumann problem
    T right_part = 0;
    T initial_distribution = 0; // Used for nonstationary and nonlinear problems

    explicit constexpr thermal_equation_data() noexcept = default;
    explicit thermal_equation_data(const Json::Value& equation) {
        energy = equation.get("energy", T{0}).template as<T>();
        right_part = equation.get("right_part", T{0}).template as<T>();
        initial_distribution = equation.get("initial_distribution", T{0}).template as<T>();
    }
};

template<std::floating_point T>
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
};

template<std::floating_point T, size_t Dimension>
struct thermal_boundaries_conditions_data final {
    std::unordered_map<std::string, thermal_boundary_condition_data<T>> conditions;

    explicit constexpr thermal_boundaries_conditions_data() noexcept = default;
    explicit thermal_boundaries_conditions_data(const Json::Value& boundaries) {
        if constexpr (Dimension == 1)
            check_required_fields(boundaries, { "left", "right" });
        for(const std::string& name : boundaries.getMemberNames())
            conditions[name] = thermal_boundary_condition_data<T>{boundaries[name]};
    }
};

template<std::floating_point T>
struct stationary_thermal_1d_data {
    Json::Value other;
    save_data save;
    thermal_equation_data<T> equation;
    thermal_boundaries_conditions_data<T, 1> boundaries;           // required
    std::vector<segment_data<T, thermal_material_data>> materials; // required
    size_t element_order = 1;
    size_t quadrature_order = 0; // 0 means that will be used default quadrature order

    static size_t convert(const std::string& order);

    explicit stationary_thermal_1d_data(const Json::Value& value)
        : other{value.get("other", {})}
        , save{value.get("save", {})} {
        check_required_fields(value, { "boundaries", "materials" });
        
        if (value.isMember("equation"))
            equation = thermal_equation_data<T>{value["equation"]};
        boundaries = thermal_boundaries_conditions_data<T, 1>{value["boundaries"]};

        const Json::Value& segments = value["materials"];
        if (!segments.isArray() || segments.empty())
            throw std::domain_error{"Field \"materials\" must be not empty array."};
        materials.reserve(segments.size());
        for(const Json::Value& segment : segments)
            materials.emplace_back(segment);

        if (value.isMember("element_order"))
            element_order = get_order(value["element_order"]);
        if (value.isMember("quadrature_order"))
            element_order = get_order(value["quadrature_order"]);
    }

    virtual ~stationary_thermal_1d_data() noexcept = default;
};

template<std::floating_point T>
struct nonstationary_thermal_1d_data final : public stationary_thermal_1d_data<T> {
    nonstationary_data<T> nonstationary;

    explicit nonstationary_thermal_1d_data(const Json::Value& value) 
        : stationary_thermal_1d_data<T>{value} {
        check_required_fields(value, { "nonstationary" });
        nonstationary = nonstationary_data<T>{value["nonstationary"]};
    }

    ~nonstationary_thermal_1d_data() noexcept override = default;
};

}

#endif
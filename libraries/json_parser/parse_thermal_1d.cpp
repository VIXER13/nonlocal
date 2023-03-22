#include "parse_thermal_1d.hpp"

#include <iostream>

namespace nonlocal::config {

model_data::model_data(const Json::Value& model) {
    check_required_fields(model, { "local_weight", "nonlocal_radius" });
    local_weight = model["local_weight"].asDouble();
    nonlocal_radius = model["nonlocal_radius"].asDouble();
    search_radius = model.get("search_radius", nonlocal_radius).asDouble();
}

physical_data::physical_data(const Json::Value& physical) {
    check_required_fields(physical, { "conductivity" });
    conductivity = physical["conductivity"].asDouble();
    capacity = physical.get("capacity", 1.).asDouble();
    density = physical.get("density", 1.).asDouble();
}

segment_data::segment_data(const Json::Value& segment) {
    check_required_fields(segment, { "elements_count", "length", "physical" });
    elements_count = segment["elements_count"].asUInt64();
    length = segment["length"].asDouble();
    if (segment.isMember("model"))
        model = model_data{segment["model"]};
    physical = physical_data{segment["physical"]};
}

thermal::boundary_condition_t thermal_boundary_condition_1d::convert(const std::string& type) {
    static const std::unordered_map<std::string, thermal::boundary_condition_t> types = {
        {"temperature", thermal::boundary_condition_t::TEMPERATURE},
        {"flux",        thermal::boundary_condition_t::FLUX},
        {"convection",  thermal::boundary_condition_t::CONVECTION},
        {"radiation",   thermal::boundary_condition_t::RADIATION},
        {"combined",    thermal::boundary_condition_t::COMBINED}
    };
    if (const auto it = types.find(type); it != types.cend())
        return it->second;
    throw std::domain_error{"Unknown boundary condition type: " + type};
}

thermal_boundary_condition_1d::thermal_boundary_condition_1d(const Json::Value& condition) {
    check_required_fields(condition, { "kind" });
    const Json::Value& type = condition["kind"];
    kind = type.isIntegral() ? thermal::boundary_condition_t(type.asUInt() - 1) : convert(type.asString());
    switch (kind) {
    case thermal::boundary_condition_t::TEMPERATURE: {
        check_required_fields(condition, { "temperature" });
        temperature = condition["temperature"].asDouble();
    } break;
    
    case thermal::boundary_condition_t::FLUX: {
        check_required_fields(condition, { "flux" });
        flux = condition["flux"].asDouble();
    } break;

    case thermal::boundary_condition_t::CONVECTION: {
        check_required_fields(condition, { "temperature", "heat_transfer" });
        temperature = condition["temperature"].asDouble();
        heat_transfer = condition["heat_transfer"].asDouble();
    } break;

    case thermal::boundary_condition_t::RADIATION: {
        check_required_fields(condition, { "emissivity" });
        emissivity = condition["emissivity"].asDouble();
    } break;

    case thermal::boundary_condition_t::COMBINED: {
        temperature = condition.get("temperature", 0.).asDouble();
        flux = condition.get("flux", 0.).asDouble();
        heat_transfer = condition.get("heat_transfer", 0.).asDouble();
        emissivity = condition.get("emissivity", 0.).asDouble();
    } break;

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(kind))};
    }
}

thermal_equation_1d_data::thermal_equation_1d_data(const Json::Value& equation) {
    check_required_fields(equation, { "boundaries" });
    const Json::Value& boundaries = equation["boundaries"];
    check_required_fields(boundaries, { "left", "right" });
    energy = equation.get("energy", 0.).asDouble();
    right_part = equation.get("right_part", 0.).asDouble();
    left = thermal_boundary_condition_1d{boundaries["left"]};
    right = thermal_boundary_condition_1d{boundaries["right"]};
}

size_t stationary_thermal_1d_data::convert(const std::string& order) {
    static const std::unordered_map<std::string, size_t> orders = {
        {"linear",    1},
        {"quadratic", 2},
        {"qubic",     3},
        {"quartic",   4},
        {"quintic",   5}
    };
    if (const auto it = orders.find(order); it != orders.cend())
        return it->second;
    throw std::domain_error("Unknown element type: " + order);  
}

stationary_thermal_1d_data::stationary_thermal_1d_data(const Json::Value& value)
    : other{value.get("other", {})}
    , save{value.get("save", {})} {
    check_required_fields(value, { "equation", "materials" });
    if (value.isMember("element_order")) {
        const Json::Value& order = value["element_order"];
        element_order = order.isIntegral() ? order.asUInt() : convert(order.asString());
        if (!element_order || element_order > 5)
            throw std::domain_error{"Invalid element order: " + std::to_string(element_order)};
    }
    equation = thermal_equation_1d_data{value["equation"]};
    const Json::Value& materials = value["materials"];
    if (!materials.isArray() || materials.empty())
        throw std::domain_error{"Field \"materials\" must be not empty array."};
    segments.reserve(materials.size());
    for(const Json::Value& segment : materials)
        segments.emplace_back(segment);
}

nonstationary_data::nonstationary_data(const Json::Value& nonstationary) {
    check_required_fields(nonstationary, { "time_step", "steps_cont"});
    time_step = nonstationary["time_step"].asDouble();
    initial_time = nonstationary.get("initial_time", 0.).asDouble();
    initial_distribution = nonstationary.get("initial_distribution", 0.).asDouble();
    steps_cont = nonstationary["steps_cont"].asUInt64();
    save_frequency = nonstationary.get("save_frequency", 1).asUInt64();
}

nonstationary_thermal_1d_data::nonstationary_thermal_1d_data(const Json::Value& value)
    : stationary_thermal_1d_data{value} {
    check_required_fields(value, { "nonstationary" });
    nonstationary = nonstationary_data{value["nonstationary"]};
}

}
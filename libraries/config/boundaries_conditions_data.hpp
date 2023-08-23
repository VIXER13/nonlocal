#ifndef NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP
#define NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP

#include "thermal_boundary_condition_data.hpp"

namespace nonlocal::config {

template<template<class, size_t> class Condition, std::floating_point T, size_t Dimension>
struct boundaries_conditions_data final {
    std::unordered_map<std::string, Condition<T, Dimension>> conditions;

    explicit constexpr boundaries_conditions_data() noexcept = default;
    explicit boundaries_conditions_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        if constexpr (Dimension == 1)
            check_required_fields(config, { "left", "right" }, path_with_access);
        for(const auto& [name, condition] : config.items()) {
            conditions[name] = Condition<T, Dimension>{condition, path_with_access + name};
        }
    }

    operator nlohmann::json() const {
        return conditions;
    }
};

template<class T>
using thermal_boundaries_conditions_1d = boundaries_conditions_data<thermal_boundary_condition_data, T, 1>;

template<class T>
using thermal_boundaries_conditions_2d = boundaries_conditions_data<thermal_boundary_condition_data, T, 2>;

}

#endif
#ifndef NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP
#define NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

template<template<class, size_t> class Condition, std::floating_point T, size_t Dimension>
struct boundaries_conditions_data final {
    std::unordered_map<std::string, Condition<T, Dimension>> conditions;

    explicit constexpr boundaries_conditions_data() noexcept = default;
    explicit boundaries_conditions_data(const nlohmann::json& boundaries) {
        if constexpr (Dimension == 1)
            check_required_fields(boundaries, { "left", "right" });
        for(const auto& [name, condition] : boundaries.items())
            conditions[name] = Condition<T, Dimension>{condition};
    }

    operator nlohmann::json() const {
        return conditions;
    }
};

}

#endif
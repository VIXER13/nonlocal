#ifndef NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP
#define NONLOCAL_CONFIG_BOUNDARIES_CONDITIONS_DATA_HPP

#include <json/value.h>

#include <type_traits>
#include <unordered_map>
#include <string>

namespace nonlocal::config {

template<template<class, size_t> class Condition, std::floating_point T, size_t Dimension>
struct boundaries_conditions_data final {
    std::unordered_map<std::string, Condition<T, Dimension>> conditions;

    explicit constexpr boundaries_conditions_data() noexcept = default;
    explicit boundaries_conditions_data(const Json::Value& boundaries) {
        if constexpr (Dimension == 1)
            check_required_fields(boundaries, { "left", "right" });
        for(const std::string& name : boundaries.getMemberNames())
            conditions[name] = Condition<T, Dimension>{boundaries[name]};
    }

    Json::Value to_json() const {
        Json::Value result;
        for(const auto& [name, condition] : conditions)
            result[name] = condition.to_json();
        return result;
    }
};

}

#endif
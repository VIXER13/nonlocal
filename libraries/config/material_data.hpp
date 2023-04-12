#ifndef NONLOCAL_CONFIG_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_MATERIAL_DATA_HPP

#include "model_data.hpp"

namespace nonlocal::config {

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct material_data final {
    Physics<T, Dimension> physical; // required
    model_data<T, Dimension> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const Json::Value& material) {
        check_required_fields(material, { "physical" });
        physical = Physics<T, Dimension>{material["physical"]};
        if (material.isMember("model"))
            model = model_data<T, Dimension>{material["model"]};
    }

    Json::Value to_json() const {
        Json::Value result;
        result["physical"] = physical.to_json();
        result["model"] = model.to_json();
        return result;
    }
};

template<template<class, size_t> class Physics, std::floating_point T>
struct material_data<Physics, T, 1> final {
    size_t elements_count = 100u; // required
    T length = T{1};              // required
    Physics<T, 1> physical;       // required
    model_data<T, 1> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const Json::Value& material) {
        check_required_fields(material, { "elements_count", "length", "physical" });
        elements_count = material["elements_count"].asUInt64();
        length = material["length"].template as<T>();
        physical = Physics<T, 1>{material["physical"]};
        if (material.isMember("model"))
            model = model_data<T, 1>{material["model"]};
    }

    Json::Value to_json() const {
        Json::Value result;
        result["elements_count"] = elements_count;
        result["length"] = length;
        result["physical"] = physical.to_json();
        result["model"] = model.to_json();
        return result;
    }
};

}

#endif
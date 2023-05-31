#ifndef NONLOCAL_CONFIG_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_MATERIAL_DATA_HPP

#include "model_data.hpp"

namespace nonlocal::config {

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct material_data final {
    Physics<T, Dimension> physical; // required
    model_data<T, Dimension> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const nlohmann::json& material) {
        check_required_fields(material, { "physical" });
        physical = Physics<T, Dimension>{material["physical"]};
        if (material.contains("model"))
            model = model_data<T, Dimension>{material["model"]};
    }

    operator nlohmann::json() const {
        return {
            {"physical", physical},
            {"model", model}
        };
    }
};

template<template<class, size_t> class Physics, std::floating_point T>
struct material_data<Physics, T, 1> final {
    size_t elements_count = 100u; // required
    T length = T{1};              // required
    Physics<T, 1> physical;       // required
    model_data<T, 1> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const nlohmann::json& material) {
        check_required_fields(material, { "elements_count", "length", "physical" });
        elements_count = material["elements_count"].get<size_t>();
        length = material["length"].get<T>();
        physical = Physics<T, 1>{material["physical"]};
        if (material.contains("model"))
            model = model_data<T, 1>{material["model"]};
    }

    operator nlohmann::json() const {
        return {
            {"elements_count", elements_count},
            {"length", length},
            {"physical", physical},
            {"model", model}
        };
    }
};

}

#endif
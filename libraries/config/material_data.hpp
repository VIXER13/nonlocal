#ifndef NONLOCAL_CONFIG_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_MATERIAL_DATA_HPP

#include "model_data.hpp"

#include <iostream>

namespace nonlocal::config {

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct material_data final {
    Physics<T, Dimension> physical; // required
    model_data<T, Dimension> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, { "physical" }, path_with_access);
        physical = Physics<T, Dimension>{config["physical"], path_with_access + "physical"};
        if (config.contains("model"))
            model = model_data<T, Dimension>{config["model"], path_with_access + "model"};
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
    explicit material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, { "elements_count", "length", "physical" }, path_with_access);
        check_optional_fields(config, {"model"}, path_with_access);
        elements_count = config["elements_count"].get<size_t>();
        length = config["length"].get<T>();
        physical = Physics<T, 1>{config["physical"], path_with_access + "physical"};
        if (config.contains("model"))
            model = model_data<T, 1>{config["model"], path_with_access + "model"};
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
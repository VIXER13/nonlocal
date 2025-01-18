#pragma once

#include "model_data.hpp"

#include <iostream>

namespace nonlocal::config {

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
std::string full_model_name() {
    if constexpr (Physics<T, Dimension>::Prefix.empty())
        return "model";
    using namespace std::string_literals;
    return Physics<T, Dimension>::Prefix.data() + "_model"s;
}

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct material_data final {
    Physics<T, Dimension> physical; // required
    model_data<T, Dimension> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, { "physical" }, path_with_access);
        physical = Physics<T, Dimension>{config["physical"], path_with_access + "physical"};
        const std::string full_name = full_model_name<Physics, T, Dimension>();
        if (const std::string model_name = config.contains(full_name) ? full_name : "model"; config.contains(model_name))
            model = model_data<T, Dimension>{config[model_name], path_with_access + model_name};
    }

    operator nlohmann::json() const {
        return {
            {"physical", physical},
            {full_model_name<Physics, T, Dimension>(), model}
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
        elements_count = config["elements_count"].get<size_t>();
        length = config["length"].get<T>();
        physical = Physics<T, 1>{config["physical"], path_with_access + "physical"};
        const std::string full_name = full_model_name<Physics,T, 1>();
        if (const std::string model_name = config.contains(full_name) ? full_name : "model"; config.contains(model_name))
            model = model_data<T, 1>{config[model_name], path_with_access + model_name};
    }

    operator nlohmann::json() const {
        using namespace std::string_literals;
        return {
            {"elements_count", elements_count},
            {"length", length},
            {"physical", physical},
            {full_model_name<Physics, T, 1>(), model}
        };
    }
};

}
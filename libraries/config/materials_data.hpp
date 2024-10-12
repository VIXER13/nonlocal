#ifndef NONLOCAL_CONFIG_THERMAL_MATERIALS_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_MATERIALS_DATA_HPP

#include "material_data.hpp"

namespace nonlocal::config {

class _materials_data final {
    template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
    friend struct materials_data;

    constexpr explicit _materials_data() noexcept = default;
    [[noreturn]] static std::string throw_error(const std::string& field, const std::string& type);
};

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct materials_data final {
    using materials_t = std::conditional_t<
        Dimension == 1,
        std::vector<material_data<Physics, T, 1>>,
        std::unordered_map<std::string, material_data<Physics, T, Dimension>>
    >;

    materials_t materials;

    explicit materials_data(const nlohmann::json& config, const std::string& path = {}) {
        using namespace std::literals;
        if constexpr (Dimension == 1) {
            if (!config.is_array() || config.empty())
                _materials_data::throw_error(path, "array");
            materials.reserve(config.size());
            for(const size_t i : std::ranges::iota_view{0u, materials.capacity()})
                materials.emplace_back(config[i], append_access_sign(path, i));
        } else {
            if (!config.is_object() || config.empty())
                _materials_data::throw_error(path, "object");
            for(const auto& [name, material] : config.items())
                materials.emplace(name, material_data<Physics, T, Dimension>{material, append_access_sign(path) + name});
        }
    }
};

}

#endif
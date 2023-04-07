#ifndef NONLOCAL_GENERAL_CONFIG_DATA_HPP
#define NONLOCAL_GENERAL_CONFIG_DATA_HPP

#include "config_utils.hpp"

#include "nonlocal_constants.hpp"

#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <optional>

namespace nonlocal::config {

class save_data final {
    std::filesystem::path _folder = ".";
    std::unordered_map<std::string, std::string> _names;
    std::optional<std::streamsize> _precision;

public:
    explicit save_data() = default;
    explicit save_data(const Json::Value& save);

    std::optional<std::streamsize> precision() const noexcept;
    const std::filesystem::path& folder() const noexcept;

    bool contains(const std::string& key) const;
    std::string get_name(const std::string& key, const std::optional<std::string>& default_name = std::nullopt) const;

    std::filesystem::path make_path(const std::string& name, const std::string& extension) const;
    std::filesystem::path path(const std::string& key, const std::string& extension,
                               const std::optional<std::string>& default_name = std::nullopt) const;

    Json::Value to_json() const;
};

template<std::floating_point T>
struct nonstationary_data final {
    T time_step = T{0.01};      // required
    T initial_time = T{0};
    uint64_t steps_cont = 100u; // required
    uint64_t save_frequency = 1u;

    explicit constexpr nonstationary_data() noexcept = default;
    explicit nonstationary_data(const Json::Value& nonstationary);

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension>
class model_data final {
    using radius_t = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;

    static radius_t read_radius(const Json::Value& arr, const std::string& field);

public:
    T local_weight = T{1};             // required
    radius_t nonlocal_radius = {T{0}}; // required
    radius_t search_radius = {T{0}};   // if skipped sets equal nonlocal_radius

    explicit constexpr model_data() noexcept = default;
    explicit model_data(const Json::Value& model);

    Json::Value to_json() const;
};

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
struct material_data final {
    Physics<T, Dimension> physical; // required
    model_data<T, Dimension> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const Json::Value& material);

    Json::Value to_json() const;
};

template<template<class, size_t> class Physics, std::floating_point T>
struct material_data<Physics, T, 1> final {
    size_t elements_count = 100u; // required
    T length = T{1};              // required
    Physics<T, 1> physical;       // required
    model_data<T, 1> model;

    explicit constexpr material_data() noexcept = default;
    explicit material_data(const Json::Value& material);

    Json::Value to_json() const;
};

template<std::floating_point T, size_t Dimension>
struct mesh_data final {
    std::filesystem::path mesh;

    explicit mesh_data() = default;
    explicit mesh_data(const Json::Value& value);

    Json::Value to_json() const;
};

template<std::floating_point T>
struct mesh_data<T, 1> final {
    size_t element_order = 1;
    size_t quadrature_order = 1;

    explicit constexpr mesh_data() noexcept = default;
    explicit mesh_data(const Json::Value& value);

    Json::Value to_json() const;
};

thermal::boundary_condition_t get_thermal_condition(const Json::Value& kind);
const std::string& get_thermal_condition(const thermal::boundary_condition_t kind);

size_t get_order(const Json::Value& order);
const std::string& get_order(const size_t order);

material_t get_material(const Json::Value& material);
const std::string& get_material(const material_t material);

}

#include "general_config_data_impl.hpp"

#endif
#pragma once

#include "config_utils.hpp"

#include <metamath/functions/norm.hpp>
#include <mesh/mesh_2d/mesh_2d.hpp>
#include <solvers/base/equation_parameters.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>

namespace nonlocal::config {

enum class influence_t : uint8_t {
    Custom,
    Constant,
    Polynomial,
    Exponential,
    Polynomial_With_Angle
};

enum class distance_t : uint8_t {
    Custom,
    Lp,
    Ellipse_With_Rotation
};

NLOHMANN_JSON_SERIALIZE_ENUM(influence_t, {
    {influence_t::Custom, nullptr},
    {influence_t::Constant, "constant"},
    {influence_t::Polynomial, "polynomial"},
    {influence_t::Exponential, "exponential"}
})

NLOHMANN_JSON_SERIALIZE_ENUM(distance_t, {
    {distance_t::Custom, nullptr},
    {distance_t::Lp, "lp"},
    {distance_t::Ellipse_With_Rotation, "ellipse_with_rotation"},
})

std::string get_model_field(const nlohmann::json& config, const std::string& path_with_access, const std::string& prefix);

class _read_model final {
    template<std::floating_point T, size_t Dimension>
    static std::array<T, Dimension> read_nonlocal_radii(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static mesh::distance_function<T> read_distance_2d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static std::function<T(T, T)> read_influence_1d(const nlohmann::json& config, const std::string& path, const T radius);

    template<std::floating_point T>
    friend std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> read_influence_2d(
        const nlohmann::json& config, const std::string& path, const std::array<T, 2>& radius);

    template<size_t Dimension, std::floating_point T>
    static bool check_parameters(const T local_weight, const std::array<T, Dimension>& radii) noexcept;

    template<size_t Dimension, std::floating_point T>
    static void fix_parameters(T& local_weight, std::array<T, Dimension>& radii) noexcept;

    explicit _read_model() noexcept = default;

public:
    template<std::floating_point T>
    friend model_parameters<1u, T> read_model_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend model_parameters<2u, T> read_model_2d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend mesh::influences<T> read_influences(const nlohmann::json& config, const std::string& path, const std::string& prefix);
};

template<std::floating_point T, size_t Dimension>
std::array<T, Dimension> _read_model::read_nonlocal_radii(const nlohmann::json& config, const std::string& path) {
    static_assert(Dimension > 0u, "The Dimension must be greater than 0.");
    std::array<T, Dimension> result;
    if (config.is_number())
        result.fill(config.get<T>());
    else if (config.is_array() && config.size() == Dimension)
        for(const size_t i : std::ranges::iota_view{0u, Dimension})
            result[i] = config[i].get<T>();
    else
        throw std::domain_error{"Field \"" + path + "\" must be an array with length " + std::to_string(Dimension)};
    return result;
}

template<std::floating_point T>
mesh::distance_function<T> read_distance_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_optional_fields(config, { "distance", "n" }, path_with_access);
    mesh::size_t_or<T> n = 2zu;
    if (config.contains("n")) {
        n = config["n"].is_number_unsigned() ? mesh::size_t_or<T>{config["n"].get<size_t>()} :
                                               mesh::size_t_or<T>{config["n"].get<T>()};
        if (mesh::get_value(n) <= T{0})
            throw std::domain_error{"Parameter \"" + path_with_access + "n\" shall be greater than 0."};
    }
    if (config.contains("distance")) {
        if (const auto distance = config["distance"].get<distance_t>(); distance == distance_t::Ellipse_With_Rotation)
            return mesh::powered_distance_with_rotation<T>{};
        else if (distance != distance_t::Lp)
            throw std::domain_error{"Unknown distance function type: " + path};
    }
    return mesh::powered_distance<T>{n};
}

template<std::floating_point T>
std::function<T(T, T)> read_influence_1d(const nlohmann::json& config, const std::string& path, const T radius) {
    using namespace nonlocal::solver_1d::influence;
    check_optional_fields(config, { "influence", "p", "q" }, path);
    if (const auto influence = config["influence"].get<influence_t>(); influence == influence_t::Constant)
        return constant_1d<T>{radius};
    else if (influence == influence_t::Exponential)
        return normal_distribution_1d<T>{radius};
    return polynomial_1d<T, 1, 1>{radius};
}

template<std::floating_point T>
std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> read_influence_2d(
    const nlohmann::json& config, const std::string& path, const std::array<T, 2>& radius) {
    using namespace nonlocal::solver_2d::influence;
    check_optional_fields(config, { "influence", "p", "q" }, path);
    const auto distance = read_distance_2d<T>(config, path);
    if (const auto influence = config["influence"].get<influence_t>(); influence == influence_t::Constant)
        return constant_2d<T>{{distance, radius}};
    else if (influence == influence_t::Exponential)
        return exponential_2d<T>{{distance, radius}, 2zu, 0.5};
    else if (influence != influence_t::Polynomial)
        throw std::domain_error{"Unknown influence function: " + path};
    return polynomial_2d<T>{{distance, radius}, 2zu, 1zu};
}

template<size_t Dimension, std::floating_point T>
bool _read_model::check_parameters(const T local_weight, const std::array<T, Dimension>& radii) noexcept {
    return (local_weight >= Nonlocal_Threshold<T> && local_weight <= T{1}) || // local problem ignores radius
           (local_weight > T{0} && local_weight < Nonlocal_Threshold<T> &&    // nonlocal problem doesn't ignore radius
            std::all_of(radii.begin(), radii.end(), [](const T radius) noexcept { return radius > T{0}; }));
}

template<size_t Dimension, std::floating_point T>
void _read_model::fix_parameters(T& local_weight, std::array<T, Dimension>& radii) noexcept {
    if (theory_type(local_weight) == theory_t::LOCAL)
        radii.fill(T{0});
    if (std::any_of(radii.begin(), radii.end(), [](const T radius) noexcept { return radius < std::numeric_limits<T>::epsilon(); }))
        local_weight = T{1};
}

template<std::floating_point T>
model_parameters<1u, T> read_model_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, { "local_weight", "nonlocal_radius" }, path_with_access);
    auto nonlocal_radius = _read_model::read_nonlocal_radii<T, 1u>(config["nonlocal_radius"], path_with_access + "nonlocal_radius");
    T local_weight = config["local_weight"].get<T>();
    if (!_read_model::check_parameters<1u>(local_weight, nonlocal_radius))
        throw std::domain_error{"Error in model parameters \"" + path + "\". "
                                "local_weight shall be in the interval (0, 1] and nonlocal_radius > 0."};
    _read_model::fix_parameters<1u>(local_weight, nonlocal_radius);
    return {
        .influence = read_influence_1d(config, path, nonlocal_radius.front()),
        .local_weight = local_weight
    };
}

template<std::floating_point T>
model_parameters<2u, T> read_model_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, { "local_weight", "nonlocal_radius" }, path_with_access);
    auto nonlocal_radius = _read_model::read_nonlocal_radii<T, 2u>(config["nonlocal_radius"], path_with_access + "nonlocal_radius");
    T local_weight = config["local_weight"].get<T>();
    if (!_read_model::check_parameters<2u>(local_weight, nonlocal_radius))
        throw std::domain_error{"Error in model parameters \"" + path + "\". "
                                "local_weight shall be in the interval (0, 1] and nonlocal_radius > 0."};
    _read_model::fix_parameters<2u>(local_weight, nonlocal_radius);
    return {
        .influence = read_influence_2d(config, path, nonlocal_radius),
        .local_weight = local_weight
    };
}

template<std::floating_point T>
mesh::influences<T> read_influences(const nlohmann::json& config, const std::string& path, const std::string& prefix) {
    mesh::influences<T> influences;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_material = append_access_sign(path) + name;
        if (const std::string model_field = get_model_field(material, path_with_material, prefix); !model_field.empty()) {
            const std::string model_path = append_access_sign(path_with_material) + model_field;
            const std::string path_with_access = append_access_sign(model_path);
            const nlohmann::json& config_model = material[model_field];
            const std::string field = config_model.contains("search_radius") ? "search_radius" :
                                      config_model.contains("nonlocal_radius") ? "nonlocal_radius" : "";
            influences[name] = {
                .distance = read_distance_2d<T>(config_model, model_path),
                .radius = field.empty() ? std::array<T, 2>{} : _read_model::read_nonlocal_radii<T, 2>(config_model[field], path_with_access + field)
            };
        }
    }
    return influences;
}

}
#pragma once

#include "config_utils.hpp"
#include "read_influence.hpp"

#include <constants/nonlocal_constants.hpp>
#include <solvers/base/equation_parameters.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>

namespace nonlocal::config {

std::string get_model_field(const nlohmann::json& config, const std::string& path_with_access, const std::string& prefix);

class _read_model final {
    template<std::floating_point T, size_t Dimension>
    static std::array<T, Dimension> read_nonlocal_radii(const nlohmann::json& config, const std::string& path);

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
    friend std::unordered_map<std::string, T> read_search_radii(const nlohmann::json& config, const std::string& path, const std::string& prefix);
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
        .influence = config.contains("influence") ? read_influence_1d(config["influence"], path_with_access + "influence", nonlocal_radius.front()) :
                                                    solver_1d::influence::polynomial_1d<T, 1u, 1u>{nonlocal_radius.front()},
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
        .influence = config.contains("influence") ? read_influence_2d(config["influence"], path_with_access + "influence", nonlocal_radius) :
                                                    solver_2d::influence::polynomial_2d<T, 2u, 1u>{nonlocal_radius},
        .local_weight = local_weight
    };
}

template<std::floating_point T>
std::unordered_map<std::string, T> read_search_radii(const nlohmann::json& config, const std::string& path, const std::string& prefix) {
    std::unordered_map<std::string, T> radii;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_material = append_access_sign(path) + name;
        if (const std::string model_field = get_model_field(material, path_with_material, prefix); !model_field.empty()) {
            const std::string path_with_access = append_access_sign(append_access_sign(path_with_material) + model_field);
            const nlohmann::json& config_model = material[model_field];
            static constexpr auto max = [](const std::array<T, 2>& radii) noexcept { return std::max(radii[0], radii[1]); };
            if (config_model.contains("search_radius"))
                radii[name] = max(_read_model::read_nonlocal_radii<T, 2u>(config_model["search_radius"], path_with_access + "search_radius"));
            else if (config_model.contains("nonlocal_radius"))
                radii[name] = max(_read_model::read_nonlocal_radii<T, 2u>(config_model["nonlocal_radius"], path_with_access + "nonlocal_radius"));
        }
    }
    return radii;
}

}
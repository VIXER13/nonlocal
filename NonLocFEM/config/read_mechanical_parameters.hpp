#pragma once

#include "read_model.hpp"

#include <solvers/solver_2d/mechanical/mechanical_parameters_2d.hpp>

#include <bitset>
#include <optional>

namespace nonlocal::config {

// return std::nullopt in isotropic case.
// return std::bitset<2> in orthotropic case.
// throw an error if parameter specified in wrong way.
std::optional<std::bitset<2>> read_null(const nlohmann::json& config, const std::string& path);

void check_null_combinations(const std::optional<std::bitset<2>> is_null_youngs_modulus,
                             const std::optional<std::bitset<2>> is_null_poissons_ratio,
                             const std::string& path);

template<std::floating_point T>
std::array<T, 2> read_elastic_parameter(const nlohmann::json& config,
                                        const std::optional<std::bitset<2>> is_null_optional) {
    std::array<T, 2> parameter;
    if (!is_null_optional)
        parameter.fill(config.get<T>());
    else {
        const auto is_null = *is_null_optional;
        parameter[0] = is_null[0] ? config[1].get<T>() : config[0].get<T>();
        parameter[1] = is_null[1] ? config[0].get<T>() : config[1].get<T>();
    }
    return parameter;
}

template<std::floating_point T>
void calculate_elastic_parameters(mechanical::parameter_2d<T>& parameters,
                                  const std::bitset<2> is_null_youngs_modulus,
                                  const std::bitset<2> is_null_poissons_ratio) noexcept {
    // youngs_modulus[1] * poissons_ratio[0] == Ey * nuxy == Ex * nuyx == youngs_modulus[0] * poissons_ratio[1]
    auto& youngs_modulus = parameters.youngs_modulus;
    auto& poissons_ratio = parameters.poissons_ratio;
    if (is_null_youngs_modulus.count() != 0) {
        if (is_null_youngs_modulus[0])
            youngs_modulus[0] = youngs_modulus[1] * poissons_ratio[0] / poissons_ratio[1];
        else
            youngs_modulus[1] = youngs_modulus[0] * poissons_ratio[1] / poissons_ratio[0];
    }
    if (is_null_poissons_ratio.count() != 0) {
        if (is_null_poissons_ratio[0])
            poissons_ratio[0] = youngs_modulus[0] * poissons_ratio[1] / youngs_modulus[1];
        else
            poissons_ratio[1] = youngs_modulus[1] * poissons_ratio[0] / youngs_modulus[0];
    }
}

template<std::floating_point T>
bool check_elasticity_parameters(const mechanical::parameter_2d<T>& parameters) noexcept {
    for(const size_t i : std::ranges::iota_view{0u, 2u}) {
        if (std::isinf(parameters.poissons_ratio[i]) || std::isinf(parameters.youngs_modulus[i]))
            return false;
        if (parameters.poissons_ratio[i] <= -1 || parameters.poissons_ratio[i] >= 0.5 || 
            std::abs(parameters.poissons_ratio[i]) < std::numeric_limits<T>::epsilon())
            return false;
        if (parameters.youngs_modulus[i] <= 0)
            return false;
    }
    if (parameters.shear_modulus <= 0)
        return false;
    return true;
}

template<std::floating_point T>
mechanical::parameter_2d<T> read_elastic_parameters(const nlohmann::json& config, const std::string& path_with_access,
                                                    const std::optional<std::bitset<2>> is_null_youngs_modulus,
                                                    const std::optional<std::bitset<2>> is_null_poissons_ratio) {
    mechanical::parameter_2d<T> parameters;
    parameters.thermal_expansion = config.value("thermal_expansion", T{0});
    parameters.youngs_modulus = read_elastic_parameter<T>(config["youngs_modulus"], is_null_youngs_modulus);
    parameters.poissons_ratio = read_elastic_parameter<T>(config["poissons_ratio"], is_null_poissons_ratio);
    if (is_null_youngs_modulus && is_null_poissons_ratio) {
        check_required_fields(config, {"shear_modulus"}, path_with_access + "shear_modulus");
        parameters.shear_modulus = config["shear_modulus"].get<T>();
        calculate_elastic_parameters(parameters, *is_null_youngs_modulus, *is_null_poissons_ratio);
    }
    return parameters;
}

template<std::floating_point T>
mechanical::parameter_2d<T> read_mechanical_coefficient_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, { "youngs_modulus", "poissons_ratio" }, path);
    check_optional_fields(config, { "thermal_expansion" }, path);
    const auto is_null_youngs_modulus = read_null(config["youngs_modulus"], path_with_access + "youngs_modulus");
    const auto is_null_poissons_ratio = read_null(config["poissons_ratio"], path_with_access + "poissons_ratio");
    check_null_combinations(is_null_youngs_modulus, is_null_poissons_ratio, path);
    const mechanical::parameter_2d<T> parameters = read_elastic_parameters<T>(config, path_with_access, is_null_youngs_modulus, is_null_poissons_ratio);
    if (!check_elasticity_parameters(parameters))
        throw std::domain_error{"Error in specifying parameters for material \"" + path + "\". "
                                "youngs_modulus must be greater than 0, "
                                "shear_modulus must be greater than 0, "
                                "shear_modulus ratio must be within the limits (-1, 0) U (0, 0.5)."};
    return parameters;
}

template<std::floating_point T>
mechanical::mechanical_parameters_2d<T> read_mechanical_parameters_2d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_object())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an object."};
    const std::string path_with_access = append_access_sign(path);
    mechanical::mechanical_parameters_2d<T> parameters;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_access_to_material = append_access_sign(path_with_access + name);
        const std::string model_field = get_model_field(material, path_with_access, "thermal");
        parameters.materials[name] = {
            .model = model_field.empty() ? model_parameters<2u, T>{} : read_model_2d<T>(material[model_field], path_with_access + model_field),
            .physical = read_mechanical_coefficient_2d<T>(material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

}
#pragma once

#include "read_model.hpp"
#include "read_coefficient.hpp"

#include <solvers/solver_2d/mechanical/mechanical_parameters_2d.hpp>

#include <bitset>
#include <optional>

namespace nonlocal::config {

class _mechanical_parameters_2d final {
    // throw an error if parameter specified in wrong way.
    static std::bitset<2> read_null(const nlohmann::json& config, const std::string& path);
    static void check_null_combinations(const std::bitset<2> is_null_young_modulus,
                                        const std::bitset<2> is_null_poissons_ratio,
                                        const std::string& path);
    template<std::floating_point T>
    static std::array<coefficient_t<T, 2>, 2> read_elastic_parameter(const nlohmann::json& config, const std::bitset<2> is_null, const std::string& path);
    template<std::floating_point T>
    static void calculate_elastic_parameters(std::array<coefficient_t<T, 2>, 2>& young_modulus, std::array<coefficient_t<T, 2>, 2>& poissons_ratio,
                                             const std::bitset<2> is_null_young_modulus, const std::bitset<2> is_null_poissons_ratio);
    template<std::floating_point T>
    static solver_2d::mechanical::isotropic_elastic_parameters<T> read_isotropic_coefficient_2d(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static solver_2d::mechanical::orthotropic_elastic_parameters<T> read_orthotropic_coefficient_2d(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static solver_2d::mechanical::anisotropic_elastic_parameters<T> read_anisotropic_coefficient_2d(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static solver_2d::mechanical::elastic_parameters_t<T> read_mechanical_coefficient_2d(const nlohmann::json& config, const std::string& path);

    explicit _mechanical_parameters_2d() noexcept = default;

public:
    template<std::floating_point T>
    friend solver_2d::mechanical::elastic_parameters<T> read_mechanical_parameters_2d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
std::array<coefficient_t<T, 2>, 2> _mechanical_parameters_2d::read_elastic_parameter(const nlohmann::json& config,
                                                                                     const std::bitset<2> is_null,
                                                                                     const std::string& path) {
    return {
        is_null[0] ? coefficient_t<T, 2>{T{0}} : read_coefficient<T, 2u>(config[0], path),
        is_null[1] ? coefficient_t<T, 2>{T{0}} : read_coefficient<T, 2u>(config[1], path)
    };
}

template<std::floating_point T>
void _mechanical_parameters_2d::calculate_elastic_parameters(std::array<coefficient_t<T, 2>, 2>& young_modulus,
                                                             std::array<coefficient_t<T, 2>, 2>& poissons_ratio,
                                                             const std::bitset<2> is_null_young_modulus,
                                                             const std::bitset<2> is_null_poissons_ratio) {
    // young_modulus[1] * poissons_ratio[0] == Ey * nuxy == Ex * nuyx == young_modulus[0] * poissons_ratio[1]
    using namespace nonlocal::utils;
    if (is_null_young_modulus.count() != 0) {
        if (is_null_young_modulus[0])
            young_modulus[0] = young_modulus[1] * poissons_ratio[0] / poissons_ratio[1];
        else
            young_modulus[1] = young_modulus[0] * poissons_ratio[1] / poissons_ratio[0];
    } else if (is_null_poissons_ratio.count() != 0) {
        if (is_null_poissons_ratio[0])
            poissons_ratio[0] = young_modulus[0] * poissons_ratio[1] / young_modulus[1];
        else
            poissons_ratio[1] = young_modulus[1] * poissons_ratio[0] / young_modulus[0];
    }
}

template<std::floating_point T>
solver_2d::mechanical::isotropic_elastic_parameters<T> _mechanical_parameters_2d::read_isotropic_coefficient_2d(
    const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    return {
        .young_modulus = read_coefficient<T, 2u>(config["young_modulus"], path_with_access + "young_modulus"),
        .poissons_ratio = read_coefficient<T, 2u>(config["poissons_ratio"], path_with_access + "poissons_ratio")
    };
}

template<std::floating_point T>
solver_2d::mechanical::orthotropic_elastic_parameters<T> _mechanical_parameters_2d::read_orthotropic_coefficient_2d(
    const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    solver_2d::mechanical::orthotropic_elastic_parameters<T> parameters = {
        .shear_modulus = read_coefficient<T, 2u>(config["shear_modulus"], path_with_access + "shear_modulus")
    };
    const auto& young_modulus = config["young_modulus"];
    const auto& poissons_ratio = config["poissons_ratio"];
    const std::string young_modulus_path = path_with_access + "young_modulus";
    const std::string poissons_ratio_path = path_with_access + "poissons_ratio";
    if (young_modulus.is_number() && poissons_ratio.is_number()) {
        parameters.young_modulus.fill(read_coefficient<T, 2u>(config["young_modulus"], young_modulus_path));
        parameters.poissons_ratio.fill(read_coefficient<T, 2u>(config["poissons_ratio"], poissons_ratio_path));
    } else if (young_modulus.is_array() && young_modulus.size() == 2 &&
               poissons_ratio.is_array() && poissons_ratio.size() == 2) {
        const auto is_null_young_modulus = read_null(config["young_modulus"], young_modulus_path);
        const auto is_null_poissons_ratio = read_null(config["poissons_ratio"], poissons_ratio_path);
        check_null_combinations(is_null_young_modulus, is_null_poissons_ratio, path);
        parameters.young_modulus = read_elastic_parameter<T>(config["young_modulus"], is_null_young_modulus, young_modulus_path);
        parameters.poissons_ratio = read_elastic_parameter<T>(config["poissons_ratio"], is_null_poissons_ratio, poissons_ratio_path);
        calculate_elastic_parameters(parameters.young_modulus, parameters.poissons_ratio,
                                     is_null_young_modulus, is_null_poissons_ratio);
    } else
        throw std::domain_error{"Wrong elastic parameters format: \"" + path + "\".\n"
                                "For orthotropic and anisotropic materials \"young_modulus\" and \"poissons_ratio\" "
                                "shall be both numbers or arrays dimension 2."};
    return parameters;
}

template<std::floating_point T>
solver_2d::mechanical::anisotropic_elastic_parameters<T> _mechanical_parameters_2d::read_anisotropic_coefficient_2d(
    const nlohmann::json& config, const std::string& path) {
    return {
        .main_parameters = read_orthotropic_coefficient_2d<T>(config, path),
        .angle = read_coefficient<T, 2u>(config["angle"], append_access_sign(path) + "angle")
    };
}

template<std::floating_point T>
solver_2d::mechanical::elastic_parameters_t<T> _mechanical_parameters_2d::read_mechanical_coefficient_2d(
    const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, { "young_modulus", "poissons_ratio" }, path);
    if (config.contains("angle")) {
        check_required_fields(config, { "shear_modulus" }, path);
        return read_anisotropic_coefficient_2d<T>(config, path);
    }
    if (config.contains("shear_modulus"))
        return read_orthotropic_coefficient_2d<T>(config, path);
    return read_isotropic_coefficient_2d<T>(config, path);
}

template<std::floating_point T>
solver_2d::mechanical::elastic_parameters<T> read_mechanical_parameters_2d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_object())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an object."};
    const std::string path_with_access = append_access_sign(path);
    solver_2d::mechanical::elastic_parameters<T> parameters;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_access_to_material = append_access_sign(path_with_access + name);
        const std::string model_field = get_model_field(material, path_with_access, "mechanical");
        parameters[name] = {
            .model = model_field.empty() ? model_parameters<2u, T>{} : read_model_2d<T>(material[model_field], path_with_access + model_field),
            .physical = _mechanical_parameters_2d::read_mechanical_coefficient_2d<T>(material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

}
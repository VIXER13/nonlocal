#pragma once

#include "read_model.hpp"
#include "read_coefficient.hpp"

#include <solvers/solver_1d/thermal/thermal_parameters_1d.hpp>
#include <solvers/solver_2d/thermal/thermal_parameters_2d.hpp>

namespace nonlocal::config {

class _read_thermal_parameters final {
    template<std::floating_point T>
    static void check_parameters(const T capacity, const T density, const T relaxation_time, const std::string& path_with_access);
    template<std::floating_point T>
    static void check_parameters(const T conductivity, const T capacity, const T density, const T relaxation_time, const std::string& path_with_access);
    // Checks only numerical parameters
    template<std::floating_point T>
    static void check_coefficient(const coefficient_t<T, 2u>& coefficient, const std::string& path_with_access);
    template<std::floating_point T>
    static void check_conductivity(const thermal::conductivity_t<T>& conductivity, const std::string& path_with_access);
    template<std::floating_point T>
    static void check_parameters(const thermal::conductivity_t<T>& conductivity, const T capacity, const T density, const T relaxation_time, const std::string& path_with_access);

    template<std::floating_point T>
    static solver_1d::thermal::parameter_1d<T> read_thermal_coefficient_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static thermal::conductivity_t<T> read_conductivity_2d(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static thermal::parameter_2d<T> read_thermal_coefficient_2d(const nlohmann::json& config, const std::string& path);

    explicit _read_thermal_parameters() noexcept = default;

public:
    template<std::floating_point T>
    friend solver_1d::thermal::parameters_1d<T> read_thermal_parameters_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend thermal::parameters_2d<T> read_thermal_parameters_2d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
void _read_thermal_parameters::check_parameters(const T capacity, const T density, const T relaxation_time, const std::string& path_with_access) {
    if (capacity <= T{0})
        throw std::domain_error{"Parameter \"" + path_with_access + "capacity\" shall be greater than 0."};
    if (density <= T{0})
        throw std::domain_error{"Parameter \"" + path_with_access + "density\" shall be greater than 0."};
    if (relaxation_time < T{0})
        throw std::domain_error{"Parameter \"" + path_with_access + "relaxation_time\" shall be greater than or equal 0."};
}

template<std::floating_point T>
void _read_thermal_parameters::check_parameters(const T conductivity, const T capacity, const T density, const T relaxation_time, const std::string& path_with_access) {
    if (conductivity <= T{0})
        throw std::domain_error{"Parameter \"" + path_with_access + "conductivity\" shall be greater than 0."};
    check_parameters(capacity, density, relaxation_time, path_with_access);
}

template<std::floating_point T>
void _read_thermal_parameters::check_coefficient(const coefficient_t<T, 2u>& coefficient, const std::string& path_with_access) {
    if (std::holds_alternative<T>(coefficient) && std::get<T>(coefficient) <= T{0})
        throw std::domain_error{"Parameter \"" + path_with_access + "conductivity\" shall be greater than 0."};
}

template<std::floating_point T>
void _read_thermal_parameters::check_conductivity(const thermal::conductivity_t<T>& conductivity, const std::string& path_with_access) {
    std::visit(metamath::visitor{
        [&](const thermal::isotropic_conductivity_t<T>& conductivity) { 
            check_coefficient(conductivity, path_with_access);
        },
        [&](const thermal::orthotropic_conductivity_t<T>& conductivity) { 
            check_coefficient(conductivity[X], path_with_access);
            check_coefficient(conductivity[Y], path_with_access);
        },
        [&](const thermal::anisotropic_conductivity_t<T>& conductivity) {
            if (std::holds_alternative<T>(conductivity[X][X]) && std::holds_alternative<T>(conductivity[X][Y]) &&
                std::holds_alternative<T>(conductivity[Y][X]) && std::holds_alternative<T>(conductivity[Y][Y])) {
                const metamath::types::square_matrix<T, 2u> matrix = {
                    std::get<T>(conductivity[X][X]), std::get<T>(conductivity[X][Y]),
                    std::get<T>(conductivity[Y][X]), std::get<T>(conductivity[Y][Y])
                };
                if (!metamath::types::is_positive(matrix))
                    throw std::domain_error{"Parameter \"" + path_with_access + "conductivity\" shall be positive matrix."};
            }
        }
    }, conductivity);
}

template<std::floating_point T>
void _read_thermal_parameters::check_parameters(
    const thermal::conductivity_t<T>& conductivity, const T capacity, const T density, const T relaxation_time, const std::string& path_with_access) {
    check_conductivity(conductivity, path_with_access);
    check_parameters(capacity, density, relaxation_time, path_with_access);
}

template<std::floating_point T>
solver_1d::thermal::parameter_1d<T> _read_thermal_parameters::read_thermal_coefficient_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"conductivity"}, path_with_access);
    check_optional_fields(config, {"capacity", "density", "relaxation_time"}, path_with_access);
    const auto [conductivity, capacity, density, relaxation_time] = std::make_tuple(
        config["conductivity"].get<T>(), config.value("capacity", T{1}), config.value("density", T{1}), config.value("relaxation_time", T{0})
    );
    check_parameters(conductivity, capacity, density, relaxation_time, path_with_access);
    return {conductivity, capacity, density, relaxation_time};
}

template<std::floating_point T>
thermal::conductivity_t<T> _read_thermal_parameters::read_conductivity_2d(const nlohmann::json& config, const std::string& path) {
    if (config.is_number())
        return read_coefficient<T, 2u>(config, path);
    if (config.is_array() && config.size() == 2)
        return thermal::orthotropic_conductivity_t<T>{ read_coefficient<T, 2u>(config[0], append_access_sign(path, 0)), 
                                                       read_coefficient<T, 2u>(config[1], append_access_sign(path, 1)) };
    if (config.is_array() && config.size() == 4)
        return thermal::anisotropic_conductivity_t<T>{
            read_coefficient<T, 2u>(config[0], append_access_sign(path, 0)),
            read_coefficient<T, 2u>(config[1], append_access_sign(path, 1)),
            read_coefficient<T, 2u>(config[2], append_access_sign(path, 2)),
            read_coefficient<T, 2u>(config[3], append_access_sign(path, 3))
        };
    throw std::domain_error{"The thermal parameter \"" + path + "\" "
                            "must be either a number in the isotropic case, "
                            "or an array of size 2 in the orthotropic case, "
                            "or an array of size 4 in the anisotropic case."};
}

template<std::floating_point T>
thermal::parameter_2d<T> _read_thermal_parameters::read_thermal_coefficient_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"conductivity"}, path_with_access);
    check_optional_fields(config, {"capacity", "density", "relaxation_time"}, path_with_access);
    const auto conductivity = read_conductivity_2d<T>(config["conductivity"], path_with_access + "conductivity");
    const auto [capacity, density, relaxation_time] = 
        std::make_tuple(config.value("capacity", T{1}), config.value("density", T{1}), config.value("relaxation_time", T{0}));
    check_parameters(conductivity, capacity, density, relaxation_time, path_with_access);
    return {conductivity, capacity, density, relaxation_time};
}

template<std::floating_point T>
solver_1d::thermal::parameters_1d<T> read_thermal_parameters_1d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_array())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an nonempty array."};

    solver_1d::thermal::parameters_1d<T> parameters(config.size());
    for(const size_t i : std::ranges::iota_view{0u, parameters.capacity()}) {
        const nlohmann::json& config_material = config[i];
        const std::string path_with_access = append_access_sign(append_access_sign(path, i));
        check_required_fields(config_material, {"physical"}, path_with_access);
        const std::string model_field = get_model_field(config_material, path_with_access, "thermal");
        using _base = _read_thermal_parameters;
        parameters[i] = {
            .model = model_field.empty() ? model_parameters<1u, T>{} : read_model_1d<T>(config_material[model_field], path_with_access + model_field),
            .physical = _base::read_thermal_coefficient_1d<T>(config_material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

template<std::floating_point T>
thermal::parameters_2d<T> read_thermal_parameters_2d(const nlohmann::json& config, const std::string& path) {
    if (!config.is_object())
        throw std::domain_error{"\"materials\" initialization requires the initializing config to be an object."};
    const std::string path_with_access = append_access_sign(path);
    thermal::parameters_2d<T> parameters;
    for(const auto& [name, material] : config.items()) {
        const std::string path_with_access_to_material = append_access_sign(path_with_access + name);
        const std::string model_field = get_model_field(material, path_with_access, "thermal");
        using _base = _read_thermal_parameters;
        parameters[name] = {
            .model = model_field.empty() ? model_parameters<2u, T>{} : read_model_2d<T>(material[model_field], path_with_access + model_field),
            .physical = _base::read_thermal_coefficient_2d<T>(material["physical"], path_with_access + "physical")
        };
    }
    return parameters;
}

}
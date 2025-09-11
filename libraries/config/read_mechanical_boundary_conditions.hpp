#pragma once

#include "config_utils.hpp"
#include "read_coefficient.hpp"

#include <logger/logger.hpp>
#include <solvers/solver_2d/mechanical/mechanical_boundary_conditions_2d.hpp>

#include <iostream>

namespace nonlocal::config {

class _mechanical_boundary_conditions final {
    template<std::floating_point T>
    static std::unique_ptr<solver_2d::mechanical::mechanical_boundary_condition_2d<T>> read_mechanical_boundary_condition_2d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static solver_2d::mechanical::mechanical_boundary_conditions_2d<T> read_mechanical_boundary_conditions_2d(const nlohmann::json& config, const std::string& path);

    explicit _mechanical_boundary_conditions() noexcept = default;

public:
    template<std::floating_point T>
    friend solver_2d::mechanical::mechanical_boundaries_conditions_2d<T> read_mechanical_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
std::unique_ptr<solver_2d::mechanical::mechanical_boundary_condition_2d<T>> 
_mechanical_boundary_conditions::read_mechanical_boundary_condition_2d(const nlohmann::json& config, const std::string& path) {
    const bool has_pressure = config.contains("pressure");
    const bool has_displacement = config.contains("displacement");
    if ((has_pressure && has_displacement) || (!has_pressure && !has_displacement))
        throw std::domain_error{"The boundary condition in \"" + path + 
                                "\" must contain only \"displacement\" or \"pressure\" field with a numerical value in it."};
    const std::string path_with_access = append_access_sign(path);
    if (has_pressure)
        return std::make_unique<solver_2d::mechanical::pressure_2d<T>>(read_coefficient<T, 2u>(config["pressure"], path_with_access + "pressure"));
    return std::make_unique<solver_2d::mechanical::displacement_2d<T>>(read_coefficient<T, 2u>(config["displacement"], path_with_access + "displacement"));
}

template<std::floating_point T>
solver_2d::mechanical::mechanical_boundary_conditions_2d<T> 
_mechanical_boundary_conditions::read_mechanical_boundary_conditions_2d(const nlohmann::json& config, const std::string& path) {
    static constexpr size_t Dimension = 2u;
    if (!config.is_array() || config.size() != Dimension)
        throw std::domain_error{"The dimension of the boundary condition \"" + path + "\" does not correspond to the dimension of the problem"};
    solver_2d::mechanical::mechanical_boundary_conditions_2d<T> conditions;
    for(const size_t i : std::ranges::iota_view{0u, Dimension}) {
        const std::string path_with_access = append_access_sign(path, i);
        if (config[i].is_null())
            logger::debug() << "The boundary condition \"" + path_with_access + "\" contain null." << std::endl;
        else
            conditions[i] = read_mechanical_boundary_condition_2d<T>(config[i], path_with_access);
    }
    return conditions;
}

template<std::floating_point T>
solver_2d::mechanical::mechanical_boundaries_conditions_2d<T> read_mechanical_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    solver_2d::mechanical::mechanical_boundaries_conditions_2d<T> boundaries_conditions;
    for(const auto& [name, conditions] : config.items())
        boundaries_conditions[name] = _mechanical_boundary_conditions::read_mechanical_boundary_conditions_2d<T>(conditions, path_with_access + name);
    return boundaries_conditions;
}

}
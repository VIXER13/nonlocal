#pragma once

#include "config_utils.hpp"

#include "mechanical/mechanical_boundary_conditions_2d.hpp"
#include <logger/logger.hpp>

namespace nonlocal::config {

enum class mechanical_boundary_condition_t : uint8_t {
    Undefined,
    Displacement,
    Pressure
};

NLOHMANN_JSON_SERIALIZE_ENUM(mechanical_boundary_condition_t, {
    {mechanical_boundary_condition_t::Undefined, nullptr},
    {mechanical_boundary_condition_t::Displacement, "displacement"},
    {mechanical_boundary_condition_t::Pressure, "pressure"}
})

template<std::floating_point T>
std::unique_ptr<mechanical::mechanical_boundary_condition_2d<T>> read_mechanical_boundary_condition_2d(const nlohmann::json& config, const std::string& path) {
    const bool has_pressure = config.contains("pressure");
    const bool has_displacement = config.contains("displacement");
    if ((has_pressure && has_displacement) || (!has_pressure && !has_displacement))
        throw std::domain_error{"The boundary condition in \"" + path + 
                                "\" must contain only \"displacement\" or \"pressure\" field with a numerical value in it."};
    if (has_pressure)
        return std::make_unique<mechanical::pressure_2d<T>>(config["pressure"].get<T>());
    return std::make_unique<mechanical::displacement_2d<T>>(config["displacement"].get<T>());
}

template<std::floating_point T>
mechanical::mechanical_boundary_conditions_2d<T> read_mechanical_boundary_conditions_2d(const nlohmann::json& config, const std::string& path) {
    static constexpr size_t Dimension = 2u;
    if (!config.is_array() || config.size() != Dimension)
        throw std::domain_error{"The dimension of the boundary condition \"" + path + "\" does not correspond to the dimension of the problem"};
    mechanical::mechanical_boundary_conditions_2d<T> conditions;
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
mechanical::mechanical_boundaries_conditions_2d<T> read_mechanical_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    mechanical::mechanical_boundaries_conditions_2d<T> boundaries_conditions;
    for(const auto& [name, conditions] : config.items())
        boundaries_conditions[name] = read_mechanical_boundary_conditions_2d<T>(conditions, path_with_access + name);
    return boundaries_conditions;
}

}
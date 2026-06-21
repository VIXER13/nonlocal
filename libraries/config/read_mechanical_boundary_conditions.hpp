#pragma once

#include "config_utils.hpp"
#include "read_coefficient.hpp"

#include <logger/logger.hpp>
#include <solvers/solver_1d/mechanical/mechanical_boundary_conditions_1d.hpp>
#include <solvers/solver_2d/mechanical/mechanical_boundary_conditions_2d.hpp>

#include <iostream>

namespace nonlocal::config {

enum class mechanical_boundary_condition_t : uint8_t {
    Undefined,
    Displacement,
    Force,
    Spring,
    Combined
};

NLOHMANN_JSON_SERIALIZE_ENUM(mechanical_boundary_condition_t, {
    {mechanical_boundary_condition_t::Undefined, nullptr},
    {mechanical_boundary_condition_t::Displacement, "displacement"},
    {mechanical_boundary_condition_t::Force,        "force"},
    {mechanical_boundary_condition_t::Spring,       "spring"},
    {mechanical_boundary_condition_t::Combined,     "combined"}
})

class _mechanical_boundary_conditions final {
    template<std::floating_point T>
    static void check_parameters(const T stiffness, const std::string& path_with_access);

    template<std::floating_point T>
    static std::unique_ptr<solver_1d::mechanical::mechanical_boundary_condition_1d<T>> read_mechanical_boundary_condition_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static std::unique_ptr<solver_2d::mechanical::mechanical_boundary_condition_2d<T>> read_mechanical_boundary_condition_2d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static solver_2d::mechanical::mechanical_boundary_conditions_2d<T> read_mechanical_boundary_conditions_2d(const nlohmann::json& config, const std::string& path);

    explicit _mechanical_boundary_conditions() noexcept = default;

public:
    template<std::floating_point T>
    friend solver_1d::mechanical::mechanical_boundaries_conditions_1d<T> read_mechanical_boundaries_conditions_1d(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend solver_2d::mechanical::mechanical_boundaries_conditions_2d<T> read_mechanical_boundaries_conditions_2d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
void _mechanical_boundary_conditions::check_parameters(const T stiffness, const std::string& path_with_access) {
    if (stiffness < T{0})
        throw std::domain_error{"\"" + path_with_access + "stiffness\" parameter shall be greather than 0."};
}

template<std::floating_point T>
std::unique_ptr<solver_1d::mechanical::mechanical_boundary_condition_1d<T>> 
_mechanical_boundary_conditions::read_mechanical_boundary_condition_1d(const nlohmann::json& config, const std::string& path) {
    using namespace solver_1d::mechanical;
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"kind"}, path_with_access);
    switch (config["kind"].get<mechanical_boundary_condition_t>()) {
    case mechanical_boundary_condition_t::Displacement:
        check_required_fields(config, { "displacement" }, path_with_access);
        return std::make_unique<displacement_1d<T>>(config["displacement"].get<T>());

    case mechanical_boundary_condition_t::Force:
        check_required_fields(config, { "force" }, path_with_access);
        return std::make_unique<normal_force_1d<T>>(config["force"].get<T>());

    case mechanical_boundary_condition_t::Spring: {
        check_required_fields(config, { "displacement", "stiffness" }, path_with_access);
        const T stiffness = config["stiffness"].get<T>();
        check_parameters(stiffness, path_with_access);
        return std::make_unique<spring_1d<T>>(stiffness, config["displacement"].get<T>());
    }

    case mechanical_boundary_condition_t::Combined: {
        if (!config.contains("stiffness"))
            check_optional_fields(config, {"force", "stiffness"}, path_with_access);
        else {
            check_required_fields(config, {"displacement"}, path_with_access);
            check_optional_fields(config, {"force"}, path_with_access);
        }
        const T stiffness = config.value("stiffness", T{0});
        check_parameters(stiffness, path_with_access);
        return std::make_unique<combined_loading_1d<T>>(
            config.value("force", T{0}),
            stiffness, config.value("displacement", T{0}));
    }

    default:
        throw std::domain_error{"Unknown boundary condition type: " + config["kind"].get<std::string>()};
    }
}

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
solver_1d::mechanical::mechanical_boundaries_conditions_1d<T> read_mechanical_boundaries_conditions_1d(const nlohmann::json& config, const std::string& path) {
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, {"left", "right"}, path_with_access);
    using _base = _mechanical_boundary_conditions;
    return {
        _base::read_mechanical_boundary_condition_1d<T>(config["left"], path_with_access + "left"),
        _base::read_mechanical_boundary_condition_1d<T>(config["right"], path_with_access + "right")
    };
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
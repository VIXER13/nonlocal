#pragma once

#include "config_utils.hpp"

#include <logger/logger.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>

#include <functional>
#include <variant>

namespace nonlocal::config {

enum class influence_family_t : uint8_t {
    Custom,
    Constant,
    Polynomial,
    Exponential,
    Polynomial_With_Angle
};

NLOHMANN_JSON_SERIALIZE_ENUM(influence_family_t, {
    {influence_family_t::Custom, nullptr},
    {influence_family_t::Constant, "constant"},
    {influence_family_t::Polynomial, "polynomial"},
    {influence_family_t::Exponential, "exponential"},
    {influence_family_t::Polynomial_With_Angle, "polynomial_with_angle"}
})

template<std::floating_point T>
std::function<T(T, T)> read_influence_1d(const nlohmann::json& config, const std::string& path, const T radius) {
    using namespace nonlocal::solver_1d::influence;
    check_required_fields(config, { "family" }, path);
    if (const auto family = config["family"].get<influence_family_t>(); family == influence_family_t::Constant)
        return constant_1d<T>{radius};
    else if (family == influence_family_t::Polynomial)
        return polynomial_1d<T, 1u, 1u>{radius};
    else if (family == influence_family_t::Exponential)
        return normal_distribution_1d<T>{radius};
    else {
        check_optional_fields(config, { "is_bounded" }, path);
        if (config.contains("is_bounded") && config["bounded"].get<bool>())
            return custom_limited_area_1d<T>{config["family"].get<std::string>(), radius};
        return custom_1d<T>{config["family"].get<std::string>()};
    }
}

template<std::floating_point T>
std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> read_influence_2d(const nlohmann::json& config, const std::string& path, const std::array<T, 2>& radius) {
    using namespace nonlocal::solver_2d::influence;
    check_required_fields(config, { "family" }, path);
    if (const auto family = config["family"].get<influence_family_t>(); family == influence_family_t::Constant)
        return constant_2d<T>{radius};
    else if (family == influence_family_t::Polynomial)
        return polynomial_2d<T, 1u, 1u>{radius};
    else if (family == influence_family_t::Exponential)
        return normal_distribution_2d<T>{radius};
    else if (family == influence_family_t::Polynomial_With_Angle)
        return polynomial_with_angle_2d<T>{radius};
    else {
        check_optional_fields(config, { "is_bounded" }, path);
        if (config.contains("is_bounded") && config["bounded"].get<bool>())
            return custom_limited_area_2d<T>{config["family"].get<std::string>(), radius};
        return custom_2d<T>{config["family"].get<std::string>()};
    }
}

}
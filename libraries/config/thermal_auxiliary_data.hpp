#pragma once

#include "config_utils.hpp"
#include "read_coefficient.hpp"

#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>

namespace nonlocal::config {

template<std::floating_point T>
struct thermal_auxiliary_data_1d final {
    T energy = 0;               // Used for Neumann problem
    T right_part = 0;
    T initial_distribution = 0; // Used for nonstationary and nonlinear problems

    explicit constexpr thermal_auxiliary_data_1d() noexcept = default;
    explicit thermal_auxiliary_data_1d(const nlohmann::json& config, const std::string& path = {}) {
        check_optional_fields(config, {"energy", "right_part", "initial_distribution"}, append_access_sign(path));
        energy = config.value("energy", T{0});
        right_part = config.value("right_part", T{0});
        initial_distribution = config.value("initial_distribution", T{0});
    }
    
    operator nlohmann::json() const {
        return {
            {"energy", energy},
            {"right_part", right_part},
            {"initial_distribution", initial_distribution}
        };
    }
};

template<std::floating_point T>
solver_2d::thermal::stationary_equation_parameters_2d<T> read_stationary_equation_parameters_2d(const nlohmann::json& config, const std::string& path = {}) {
    check_optional_fields(config, {"right_part", "initial_distribution", "energy", "tolerance", "max_iterations"}, append_access_sign(path));
    solver_2d::thermal::stationary_equation_parameters_2d<T> result;
    if (config.contains("right_part"))
        result.right_part = [right_part = read_coefficient<T, 2u>(config["right_part"], path)](const std::array<T, 2>& x) {
            return std::visit(metamath::types::visitor{
                [](const T value) { return value; },
                [&x](const spatial_dependency<T, 2>& value) { return value(x); },
                [](const auto&) { throw std::domain_error{"Unsuported right part format."}; return T{0}; }
            }, right_part);
        };
    if (config.contains("initial_distribution"))
        result.initial_distribution = [initial_distribution = read_coefficient<T, 2u>(config["initial_distribution"], path)](const std::array<T, 2>& x) {
            return std::visit(metamath::types::visitor{
                [](const T value) { return value; },
                [&x](const spatial_dependency<T, 2>& value) { return value(x); },
                [](const auto&) { throw std::domain_error{"Unsuported right part format."}; return T{0}; }
            }, initial_distribution);
        };
    result.energy = config.value("energy", result.energy);
    result.tolerance = config.value("tolerance", result.tolerance);
    result.max_iterations = config.value("max_iterations", result.max_iterations);
    return result;
}

}
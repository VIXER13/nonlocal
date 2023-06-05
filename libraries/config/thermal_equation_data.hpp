#ifndef NONLOCAL_CONFIG_THERMAL_EQUATION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_EQUATION_DATA_HPP

#include <nlohmann/json.hpp>

namespace nonlocal::config {

template<std::floating_point T>
struct thermal_equation_data final {
    T energy = 0;               // Used for Neumann problem
    T right_part = 0;
    T initial_distribution = 0; // Used for nonstationary and nonlinear problems

    explicit constexpr thermal_equation_data() noexcept = default;
    explicit thermal_equation_data(const nlohmann::json& equation) {
        energy = equation.value("energy", T{0});
        right_part = equation.value("right_part", T{0});
        initial_distribution = equation.value("initial_distribution", T{0});
    }

    operator nlohmann::json() const {
        return {
            {"energy", energy},
            {"right_part", right_part},
            {"initial_distribution", initial_distribution}
        };
    }
};

}

#endif
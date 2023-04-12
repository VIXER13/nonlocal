#ifndef NONLOCAL_CONFIG_THERMAL_EQUATION_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_EQUATION_DATA_HPP

namespace nonlocal::config {

template<std::floating_point T>
struct thermal_equation_data final {
    T energy = 0;               // Used for Neumann problem
    T right_part = 0;
    T initial_distribution = 0; // Used for nonstationary and nonlinear problems

    explicit constexpr thermal_equation_data() noexcept = default;
    explicit thermal_equation_data(const Json::Value& equation) {
        energy = equation.get("energy", T{0}).template as<T>();
        right_part = equation.get("right_part", T{0}).template as<T>();
        initial_distribution = equation.get("initial_distribution", T{0}).template as<T>();
    }

    Json::Value to_json() const {
        Json::Value result;
        result["energy"] = energy;
        result["right_part"] = right_part;
        result["initial_distribution"] = initial_distribution;
        return result;
    }
};

}

#endif
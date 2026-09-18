#pragma once

#include "config_utils.hpp"
#include "read_coefficient.hpp"

namespace nonlocal::config {

template<std::floating_point T>
struct mechanical_auxiliary_data_1d final {
    coefficient_t<T, 1u> right_part;
    coefficient_t<T, 1u> initial_distribution; // Used for nonstationary and nonlinear problems

    explicit constexpr mechanical_auxiliary_data_1d() noexcept = default;
    explicit mechanical_auxiliary_data_1d(const nlohmann::json& config, const std::string& path = {}) {
        check_optional_fields(config, {"right_part", "initial_distribution"}, append_access_sign(path));
        right_part           = read_coefficient<T, 1u>(config["right_part"], path);
        initial_distribution = read_coefficient<T, 1u>(config["initial_distribution"], path);
    }
    
    operator nlohmann::json() const {
        return {
            {"right_part", right_part},
            {"initial_distribution", initial_distribution}
        };
    }
};

}
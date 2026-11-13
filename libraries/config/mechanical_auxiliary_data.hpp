#pragma once

#include "config_utils.hpp"
#include "read_coefficient.hpp"

namespace nonlocal::config {

template<std::floating_point T>
struct mechanical_auxiliary_data_1d final {
    coefficient_t<T, 1u> right_part = T{ 0 };
    coefficient_t<T, 1u> initial_distribution = T{ 0 }; // Used for nonstationary and nonlinear problems

    explicit constexpr mechanical_auxiliary_data_1d() noexcept = default;
    explicit mechanical_auxiliary_data_1d(const nlohmann::json& config, const std::string& path = {}) {
        if (!config.is_object())
            throw std::domain_error{ "Field \"" + path + "\" must be an object." };
        const std::string path_with_access = append_access_sign(path);
        check_optional_fields(config, { "right_part", "initial_distribution" }, path_with_access);
        if (config.contains("right_part"))
            right_part = read_coefficient<T, 1u>(config["right_part"], path_with_access + "right_part");
        if (config.contains("initial_distribution"))
            initial_distribution = read_coefficient<T, 1u>(config["initial_distribution"], path_with_access + "initial_distribution");
    }

    operator nlohmann::json() const {
        return { { "right_part", right_part }, { "initial_distribution", initial_distribution } };
    }
};

} // namespace nonlocal::config
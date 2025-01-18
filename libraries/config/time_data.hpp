#pragma once

#include "config_utils.hpp"

namespace nonlocal::config {

template<std::floating_point T>
struct time_data final {
    T time_step = T{0};       // required
    T initial_time = T{0};
    uint64_t steps_count = 0; // required
    uint64_t save_frequency = 1ull;

    explicit constexpr time_data() noexcept = default;
    explicit time_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, {"time_step", "steps_count"}, path_with_access);
        check_optional_fields(config, {"initial_time", "save_frequency"}, path_with_access);
        time_step = config["time_step"].get<T>();
        initial_time = config.value("initial_time", T{0});
        steps_count = config["steps_count"].get<uint64_t>();
        save_frequency = config.value("save_frequency", 1ull);
    }

    operator nlohmann::json() const {
        return {
            {"time_step", time_step},
            {"initial_time", initial_time},
            {"steps_count", steps_count},
            {"save_frequency", save_frequency}
        };
    }
};

}
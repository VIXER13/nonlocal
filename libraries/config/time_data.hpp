#ifndef NONLOCAL_CONFIG_TIME_DATA_HPP
#define NONLOCAL_CONFIG_TIME_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

template<std::floating_point T>
struct time_data final {
    T time_step = T{0};       // required
    T initial_time = T{0};
    uint64_t steps_count = 0; // required
    uint64_t save_frequency = 1ull;

    explicit constexpr time_data() noexcept = default;
    explicit time_data(const nlohmann::json& time) {
        check_required_fields(time, { "time_step", "steps_count" });
        time_step = time["time_step"].get<T>();
        initial_time = time.value("initial_time", T{0});
        steps_count = time["steps_count"].get<uint64_t>();
        save_frequency = time.value("save_frequency", 1ull);
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

#endif
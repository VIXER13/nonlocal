#ifndef NONLOCAL_CONFIG_TIME_DATA_HPP
#define NONLOCAL_CONFIG_TIME_DATA_HPP

namespace nonlocal::config {

template<std::floating_point T>
struct time_data final {
    T time_step = T{0};       // required
    T initial_time = T{0};
    uint64_t steps_count = 0; // required
    uint64_t save_frequency = 1u;

    explicit constexpr time_data() noexcept = default;
    explicit time_data(const Json::Value& nonstationary) {
        check_required_fields(nonstationary, { "time_step", "steps_count"});
        time_step = nonstationary["time_step"].template as<T>();
        initial_time = nonstationary.get("initial_time", T{0}).template as<T>();
        steps_count = nonstationary["steps_count"].asUInt64();
        save_frequency = nonstationary.get("save_frequency", 1u).asUInt64();
    }

    Json::Value to_json() const {
        Json::Value result;
        result["time_step"] = time_step;
        result["initial_time"] = initial_time;
        result["steps_count"] = steps_count;
        result["save_frequency"] = save_frequency;
        return result;
    }
};

}

#endif
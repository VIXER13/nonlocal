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
    explicit time_data(const Json::Value& time) {
        check_required_fields(time, { "time_step", "steps_count"});
        time_step = time["time_step"].template as<T>();
        initial_time = time.get("initial_time", T{0}).template as<T>();
        steps_count = time["steps_count"].asUInt64();
        save_frequency = time.get("save_frequency", 1u).asUInt64();
    }

    operator Json::Value() const {
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
#pragma once

#include "config_utils.hpp"

namespace nonlocal::config {

enum class frequency_units_t : uint8_t {
    Hz,
    KHz,
    MHz,
    GHz,
    Unknown
};

NLOHMANN_JSON_SERIALIZE_ENUM(frequency_units_t, {
    {frequency_units_t::Unknown, nullptr},
    {frequency_units_t::Hz,  "Hz"},
    {frequency_units_t::KHz, "KHz"},
    {frequency_units_t::MHz, "MHz"},
    {frequency_units_t::GHz, "GHz"}
})

template<std::floating_point T>
class frequency_data final {
public:
    std::vector<T> frequencies = { T{0} }; // required

private:
    void check_frequency(const T frequency) {
        if (frequency < T{0})
            throw std::domain_error{"Frequency must be greater than or equal to zero."};
    }

    void check_min_max_frequencies(const T min_frequency, const T max_frequency) {
        check_frequency(min_frequency);
        check_frequency(max_frequency);
        if (max_frequency < min_frequency)
            throw std::domain_error{"Maximal frequency must be greater than minimal frequency."};
    }

    T frequency_units_multiplier(const frequency_units_t units) {
        switch (units) {
        case frequency_units_t::Hz:  return static_cast<T>(1e0);
        case frequency_units_t::KHz: return static_cast<T>(1e3);
        case frequency_units_t::MHz: return static_cast<T>(1e6);
        case frequency_units_t::GHz: return static_cast<T>(1e9);
        case frequency_units_t::Unknown:
        default: throw std::domain_error{"Unsupported frequency units."};
        }
    }
    
public:
    explicit constexpr frequency_data() noexcept = default;
    explicit frequency_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        if (config.contains("min_value") and config.contains("max_value") and config.contains("number_of_points")) {
            check_required_fields(config, {"units", "min_value", "max_value", "number_of_points"}, path_with_access);
            const T multiplier = frequency_units_multiplier(config["units"].get<frequency_units_t>());
            const T min_frequency = config["min_value"].get<T>() * multiplier;
            const T max_frequency = config["max_value"].get<T>() * multiplier;
            check_min_max_frequencies(min_frequency, max_frequency);
            const uint64_t number_of_points = config["number_of_points"].get<uint64_t>();
            const T frequency_step = (max_frequency - min_frequency) / number_of_points;
            frequencies.resize(number_of_points);
            for(const size_t i : std::ranges::iota_view{0zu, number_of_points})
                frequencies[i] = static_cast<T>(i) * frequency_step + min_frequency;
        } else if (config.contains("range")) {
            frequencies = config["range"].get<std::vector<T>>();
            for(const T frequency : frequencies) check_frequency(frequency);
        } else {
            throw std::domain_error{"Unsupported frequency parameters in \"" + path +
                                    "\". Specify either {units, min_value, max_value, number_of_points} or {range}."};
        }
    }

    operator nlohmann::json() const {
        return {
            {"frequencies", frequencies}
        };
    }

};

}

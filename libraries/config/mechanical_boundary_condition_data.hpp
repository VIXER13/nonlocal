#ifndef MECHANICAL_BOUNDARY_CONDITION_DATA_HPP
#define MECHANICAL_BOUNDARY_CONDITION_DATA_HPP

#include "config_utils.hpp"

#include "logger.hpp"

namespace nonlocal::config {

enum class mechanical_boundary_condition_t : uint8_t {
    UNDEFINED,
    DISPLACEMENT,
    PRESSURE
};

NLOHMANN_JSON_SERIALIZE_ENUM(mechanical_boundary_condition_t, {
    {mechanical_boundary_condition_t::UNDEFINED, nullptr},
    {mechanical_boundary_condition_t::DISPLACEMENT, "displacement"},
    {mechanical_boundary_condition_t::PRESSURE, "pressure"}
})

template<std::floating_point T>
struct mechanical_boundary_condition_component_data final {
    mechanical_boundary_condition_t kind = mechanical_boundary_condition_t::UNDEFINED;
    T value = T{0};
};

template<std::floating_point T, size_t Dimension>
class mechanical_boundary_condition_data final {
    using boundary_component = mechanical_boundary_condition_component_data<T>;

    boundary_component get_condition(const nlohmann::json& config, const std::string& path) {
        if (const bool has_pressure = config.contains("pressure"); !has_pressure && !config.contains("displacement"))
            logger::get().log(logger::log_level::WARNING) << "No boundary condition type specified \"" + path + "\", so default is 0 pressure." << std::endl;
        else {
            using enum mechanical_boundary_condition_t;
            return has_pressure ? boundary_component{.kind = PRESSURE, .value = config["pressure"].get<T>()} :
                                  boundary_component{.kind = DISPLACEMENT, .value = config["displacement"].get<T>()};
        }
        return {};
    }

public:
    std::conditional_t<Dimension == 1,
        boundary_component,
        std::array<boundary_component, Dimension>
    > condition;

    explicit constexpr mechanical_boundary_condition_data() noexcept = default;
    explicit mechanical_boundary_condition_data(const nlohmann::json& config, const std::string& path = {}) {
        const nlohmann::json& conf = config.contains("mechanical") ? config["mechanical"] : config;
        const std::string& config_path = config.contains("mechanical") ? append_access_sign(path) + "mechanical" : path;
        if constexpr (Dimension == 1) 
            condition = get_condition(conf, config_path);
        else {
            if (!conf.is_array() || conf.size() != Dimension)
                throw std::domain_error{"The dimension of the boundary condition \"" + config_path + "\" does not correspond to the dimension of the problem"};
            for(const size_t i : std::ranges::iota_view{0u, Dimension}) 
                condition[i] = get_condition(conf[i], append_access_sign(config_path, i));
        }
    }
};

}

#endif
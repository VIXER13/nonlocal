#pragma once

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
struct mechanical_boundary_condition_data final {
    mechanical_boundary_condition_t kind = mechanical_boundary_condition_t::PRESSURE;
    T value = T{0};

    explicit constexpr mechanical_boundary_condition_data() noexcept = default;
    explicit mechanical_boundary_condition_data(const nlohmann::json& config, const std::string& path = {}) {
        const bool has_pressure = config.contains("pressure");
        const bool has_displacement = config.contains("displacement");
        if ((has_pressure && has_displacement) || (!has_pressure && !has_displacement))
            throw std::domain_error{"The boundary condition in \"" + path + 
                                    "\" must contain only \"displacement\" or \"pressure\" field with a numerical value in it."};
        else if (has_pressure)
            value = config["pressure"].get<T>();
        else {
            kind = mechanical_boundary_condition_t::DISPLACEMENT;
            value = config["displacement"].get<T>();
        }
    }
};

template<std::floating_point T, size_t Dimension>
struct mechanical_boundary_conditions_data final {
    std::conditional_t<Dimension == 1,
        mechanical_boundary_condition_data<T>,
        std::array<std::optional<mechanical_boundary_condition_data<T>>, Dimension>
    > conditions;

    explicit constexpr mechanical_boundary_conditions_data() noexcept = default;
    explicit mechanical_boundary_conditions_data(const nlohmann::json& config, const std::string& path = {}) {
        if constexpr (Dimension == 1) 
            conditions = mechanical_boundary_condition_data<T>{config, path};
        else {
            if (!config.is_array() || config.size() != Dimension)
                throw std::domain_error{"The dimension of the boundary condition \"" + path + "\" does not correspond to the dimension of the problem"};
            for(const size_t i : std::ranges::iota_view{0u, Dimension}) {
                const std::string path_with_access = append_access_sign(path, i);
                if (!config[i].is_null())
                    conditions[i] = mechanical_boundary_condition_data<T>{config[i], path_with_access};
                else
                    logger::debug() << "The boundary condition \"" + path_with_access + "\" contain null." << std::endl;
            }
        }
    }
};

}
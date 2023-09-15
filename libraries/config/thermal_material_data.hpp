#ifndef NONLOCAL_CONFIG_THERMAL_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_MATERIAL_DATA_HPP

#include "config_utils.hpp"

#include "nonlocal_constants.hpp"

#include <type_traits>
#include <exception>
#include <ranges>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
struct thermal_material_data;

template<std::floating_point T>
struct thermal_material_data<T, 1> final {
    T conductivity = T{1}; // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string right_part = append_access_sign(path);
        check_required_fields(config, { "conductivity" }, right_part);
        check_optional_fields(config, {"capacity", "density"}, right_part);
        conductivity = config["conductivity"].get<T>();
        capacity = config.value("capacity", T{1});
        density = config.value("density", T{1});
    }

    operator nlohmann::json() const {
        return {
            {"conductivity", conductivity},
            {"capacity", capacity},
            {"density", density}
        };
    }
};

template<std::floating_point T>
class thermal_material_data<T, 2> final {
    void read_conductivity(const nlohmann::json& conduct) {
        if (conduct.is_number())
            conductivity.front() = conduct.get<T>();
        else if (conduct.is_array() && conduct.size() == 2) {
            material = material_t::ORTHOTROPIC;
            conductivity.front() = conduct[0].get<T>();
            conductivity.back() = conduct[1].get<T>();
        } else if (conduct.is_array() && conduct.size() == 4) {
            material = material_t::ANISOTROPIC;
            for(const size_t i : std::ranges::iota_view{0u, 4u})
                conductivity[i] = conduct[i].get<T>();
        } else
            throw std::domain_error{"The thermal conductivity should be either a number in the isotropic case, "
                                    "or an array of size 2 in the orthotropic case and size 4 in the anisotropic case"};
    }

public:
    material_t material = material_t::ISOTROPIC; // not json field
    std::array<T, 4> conductivity = {T{1}};      // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string right_part = append_access_sign(path);
        check_required_fields(config, { "conductivity" }, right_part);
        check_optional_fields(config, {"capacity", "density"}, right_part);
        read_conductivity(config["conductivity"]);
        capacity = config.value("capacity", T{1});
        density = config.value("density", T{1});
    }

    operator nlohmann::json() const {
        nlohmann::json result = {
            {"capacity", capacity},
            {"density", density}
        };
        if (material == material_t::ISOTROPIC)
            result["conductivity"] = conductivity.front();
        else if (material == material_t::ORTHOTROPIC)
            result["conductivity"] = std::array{conductivity.front(), conductivity.back()};
        else
            result["conductivity"] = conductivity;
        return result;
    }
};

}

#endif
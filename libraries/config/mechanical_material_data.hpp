#ifndef NONLOCAL_CONFIG_MECHANICAL_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_MECHANICAL_MATERIAL_DATA_HPP

#include "config_utils.hpp"

#include "nonlocal_constants.hpp"

#include <type_traits>
#include <exception>
#include <ranges>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
struct mechanical_material_data;

template<std::floating_point T>
struct mechanical_material_data<T, 1> final {
    static constexpr std::string_view Prefix = "mechanical";

    T youngs_modulus = T{1}; // required

    explicit constexpr mechanical_material_data() noexcept = default;
    explicit mechanical_material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string right_part = append_access_sign(path);
        check_required_fields(config, { "youngs_modulus" }, right_part);
        youngs_modulus = config["youngs_modulus"].get<T>();
    }

    operator nlohmann::json() const {
        return {
            {"youngs_modulus", youngs_modulus}
        };
    }
};

template<std::floating_point T>
struct mechanical_material_data<T, 2> final {
    static constexpr std::string_view Prefix = "mechanical";

    T youngs_modulus = T{1}; // required
    T poissons_ratio = T{0.3}; // required
    T thermal_expansion = T{0}; // optional

    explicit constexpr mechanical_material_data() noexcept = default;
    explicit mechanical_material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string right_part = append_access_sign(path);
        check_required_fields(config, { "youngs_modulus", "poissons_ratio" }, right_part);
        youngs_modulus = config["youngs_modulus"].get<T>();
        poissons_ratio = config["poissons_ratio"].get<T>();
        thermal_expansion = config.value("thermal_expansion", T{0});
    }

    operator nlohmann::json() const {
        return {
            {"youngs_modulus", youngs_modulus},
            {"poissons_ratio", poissons_ratio},
            {"thermal_expansion", thermal_expansion}
        };
    }
};

}

#endif
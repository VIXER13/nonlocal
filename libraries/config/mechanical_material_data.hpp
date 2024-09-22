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
class mechanical_material_data<T, 2> final {
    void read_elasticity_parameters(std::array<T, 2>& parameter, const nlohmann::json& config, std::string name) {
        const auto& conf = config[name];
        if (config.is_number()) {
            parameter[0] = conf.get<T>();
            parameter[1] = conf.get<T>();
        }
        else if (conf.is_array() && conf.size() == 2) {
            material = material_t::ORTHOTROPIC;
            parameter[0] = conf[0].get<T>();
            parameter[1] = conf[1].get<T>();
        } else if (conf.is_array() && conf.size() == 4) {
            throw std::domain_error{"Anisotropic case is not implemented"};
        } else
            throw std::domain_error{"The " + name + " should be either a number in the isotropic case, "
                                    "or an array of size 2 in the orthotropic"};
    }

public:
    static constexpr std::string_view Prefix = "mechanical";
    material_t material = material_t::ISOTROPIC;        // not json field
    std::array<T, 2> youngs_modulus = {T{1},   T{1}};   // required
    std::array<T, 2> poissons_ratio = {T{0.3}, T{0.3}}; // required
    T shear_modulus = T{1};
    T thermal_expansion = T{0}; // optional

    explicit constexpr mechanical_material_data() noexcept = default;
    explicit mechanical_material_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string right_part = append_access_sign(path);
        check_required_fields(config, { "youngs_modulus", "poissons_ratio"}, right_part);
        check_optional_fields(config, { "shear_modulus" }, right_part); 
        read_elasticity_parameters(youngs_modulus, config, "youngs_modulus");
        read_elasticity_parameters(poissons_ratio, config, "poissons_ratio");
        shear_modulus = config.value("shear_modulus", T{1});
        thermal_expansion = config.value("thermal_expansion", T{0});
    }

    operator nlohmann::json() const {
        return {
            {"youngs_modulus", youngs_modulus},
            {"poissons_ratio", poissons_ratio},
            {"shear_modulus", shear_modulus},
            {"thermal_expansion", thermal_expansion}
        };
    }
};

}

#endif
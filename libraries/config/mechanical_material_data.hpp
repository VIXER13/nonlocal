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
    void read_elasticity_parameters(std::array<T, 2>& parameter, const nlohmann::json& config, std::string name, std::array<size_t, 2>& null) {
        if (const auto& conf = config[name]; conf.is_number()) {
            parameter.fill(conf.get<T>());
            null[1] = 1;
        }
        else if (conf.is_array() && conf.size() == 2) {
            material = material_t::ANISOTROPIC;
            if(!conf[0].is_null() && !conf[1].is_null()) {
                parameter[0] = conf[0].get<T>();
                parameter[1] = conf[1].get<T>();
            } else if (conf[0].is_null() && !conf[1].is_null()) {
                parameter[0] = conf[1].get<T>();
                parameter[1] = conf[1].get<T>();
                null[0] = 1;
            } else if (!conf[0].is_null() && conf[1].is_null()) {
                parameter[0] = conf[0].get<T>();
                parameter[1] = conf[0].get<T>();
                null[1] = 1;
            } else {
                null[0] = 1; null[1] = 1;
            }
        } else
            throw std::domain_error{"The " + name + " should be either a number in the isotropic case, "
                                    "or an array of size 2 in the anisotropic"};
    }
    void calculate_elasticity_parameters(const std::array<size_t, 2>& null_E, const std::array<size_t, 2>& null_nu) {
        // youngs_modulus[1] * poissons_ratio[0] = Ey * nuxy = Ex * nuyx = youngs_modulus[0] * poissons_ratio[1]
        if (null_E[0] == 0 && null_E[1] == 1) {
            youngs_modulus[1] = youngs_modulus[0] * poissons_ratio[1] / poissons_ratio[0];
            return;
        }
        if (null_E[0] == 1 && null_E[1] == 0) {
            youngs_modulus[0] = youngs_modulus[1] * poissons_ratio[0] / poissons_ratio[1];
            return;
        }
        if (null_nu[0] == 0 && null_nu[1] == 1) {
            poissons_ratio[1] = youngs_modulus[1] * poissons_ratio[0] / youngs_modulus[0];
            return;
        }
        if (null_nu[0] == 1 && null_nu[1] == 0) {
            poissons_ratio[0] = youngs_modulus[0] * poissons_ratio[1] / youngs_modulus[1];
            return;
        }
    }  
    bool check_elasticity_parameters(const std::array<size_t, 2>& null_E, const std::array<size_t, 2>& null_nu) {
        const size_t null = null_nu[0] + null_nu[1] + null_E[0] + null_E[1];
        if (null == 1) return false;
        return true;
    }
    bool elasticity_parameters_correctness() {
        if (poissons_ratio[0] <= 0 || poissons_ratio[0] >= 0.5)
            return true;
        if (poissons_ratio[1] <= 0 || poissons_ratio[1] >= 0.5)
            return true;
        if (youngs_modulus[0] <= 0 || youngs_modulus[1] <= 0)
            return true;
        if (shear_modulus <= 0)
            return true;
        return false;
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
        std::array<size_t, 2> null_E {0, 0}, null_nu {0, 0};
        read_elasticity_parameters(youngs_modulus, config, "youngs_modulus", null_E);
        read_elasticity_parameters(poissons_ratio, config, "poissons_ratio", null_nu);
        if (check_elasticity_parameters(null_E, null_nu))
            throw std::domain_error{"Input format error. Choose three of the input parameters : Ex, Ey, nuxy, nuyx. The fourth will be calculated automatically"};
        if (elasticity_parameters_correctness())
            throw std::domain_error{"Input format error. All elasticity parameters must be positive and nonzero. Moreover : 0 < poissons_ratio < 0.5"};
        calculate_elasticity_parameters(null_E, null_nu);
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
#ifndef NONLOCAL_CONFIG_MECHANICAL_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_MECHANICAL_MATERIAL_DATA_HPP

#include "config_utils.hpp"

#include "nonlocal_constants.hpp"

#include <type_traits>
#include <exception>
#include <ranges>
#include <bitset>

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

    void read_elasticity_parameters(const nlohmann::json& config) {
        auto fill_null = [&](std::string name) {
            std::bitset<2> null{"00"};
            if (const auto& conf = config[name]; conf.is_number()) {
                null.flip(1);
            }
            else if (conf.is_array() && conf.size() == 2) {
                if (conf[0].is_null()) null.flip(0);
                if (conf[1].is_null()) null.flip(1);
            } else 
                throw std::domain_error{"The " + name + " should be either a number in the isotropic case, " +
                                        "or an array of size 2 in the anisotropic"};
            return null;
        };
        const std::bitset<2> null_E = fill_null("youngs_modulus"), null_nu = fill_null("poissons_ratio");
        if (!(null_E.count() + null_nu.count() == 1))
            throw std::domain_error{"Input format error. Choose three of the input parameters : Ex, Ey, nuxy, nuyx. The fourth will be calculated automatically"};
        auto read_params = [&](std::string name, const std::bitset<2>& null) {
            std::array<T, 2> parameter {0, 0};
            if (const auto& conf = config[name]; conf.is_number()) {
                parameter.fill(conf.get<T>());
            }
            else if (conf.is_array() && conf.size() == 2) {
                material = material_t::ANISOTROPIC;
                parameter[0] = null[0] ? conf[1].get<T>() : conf[0].get<T>();
                parameter[1] = null[1] ? conf[0].get<T>() : conf[1].get<T>();
            } else
                throw std::domain_error{"The " + name + " should be either a number in the isotropic case, "
                                        "or an array of size 2 in the anisotropic"};
            return parameter;
        };
        youngs_modulus = read_params("youngs_modulus", null_E);
        poissons_ratio = read_params("poissons_ratio", null_nu);
        if (!check_elasticity_parameters())
            throw std::domain_error{"Input format error. All elasticity parameters must be nonzero. Moreover : poissons_ratio in [-1, 0) U (0, 0.5]"};
        calculate_elasticity_parameters(null_E, null_nu);
    }

    void calculate_elasticity_parameters(const std::bitset<2>& null_E, const std::bitset<2>& null_nu) {
        // youngs_modulus[1] * poissons_ratio[0] = Ey * nuxy = Ex * nuyx = youngs_modulus[0] * poissons_ratio[1]
        if (!null_E[0] && null_E[1]) {
            youngs_modulus[1] = youngs_modulus[0] * poissons_ratio[1] / poissons_ratio[0];
            return;
        }
        if (null_E[0] && !null_E[1]) {
            youngs_modulus[0] = youngs_modulus[1] * poissons_ratio[0] / poissons_ratio[1];
            return;
        }
        if (!null_nu[0] && null_nu[1]) {
            poissons_ratio[1] = youngs_modulus[1] * poissons_ratio[0] / youngs_modulus[0];
            return;
        }
        if (null_nu[0] && !null_nu[1]) {
            poissons_ratio[0] = youngs_modulus[0] * poissons_ratio[1] / youngs_modulus[1];
            return;
        }
    }  

    bool check_elasticity_parameters() const noexcept {
        if (poissons_ratio[0] <= -1 || poissons_ratio[0] >= 0.5 || std::abs(poissons_ratio[0]) < 1e-16)
            return false;
        if (poissons_ratio[1] <= -1 || poissons_ratio[1] >= 0.5 || std::abs(poissons_ratio[1]) < 1e-16)
            return false;
        if (youngs_modulus[0] <= 0 || youngs_modulus[1] <= 0)
            return false;
        if (shear_modulus <= 0)
            return false;
        return true;
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
        read_elasticity_parameters(config);
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
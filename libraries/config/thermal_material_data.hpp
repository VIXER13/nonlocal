#ifndef NONLOCAL_CONFIG_THERMAL_MATERIAL_DATA_HPP
#define NONLOCAL_CONFIG_THERMAL_MATERIAL_DATA_HPP

#include <json/value.h>

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
    explicit thermal_material_data(const Json::Value& physical) {
        check_required_fields(physical, { "conductivity" });
        conductivity = physical["conductivity"].template as<T>();
        capacity = physical.get("capacity", T{1}).template as<T>();
        density = physical.get("density", T{1}).template as<T>();
    }

    Json::Value to_json() const {
        Json::Value result;
        result["conductivity"] = conductivity;
        result["capacity"] = capacity;
        result["density"] = density;
        return result;
    }
};

template<std::floating_point T>
class thermal_material_data<T, 2> final {
    void read_conductivity(const Json::Value& conduct) {
        if (conduct.isNumeric())
            conductivity.front() = conduct.template as<T>();
        else if (conduct.isArray() && conduct.size() == 2) {
            material == material_t::ORTHOTROPIC;
            conductivity.front() = conduct[0].template as<T>();
            conductivity.back() = conduct[1].template as<T>();
        } else if (conduct.isArray() && conduct.size() == 4) {
            material == material_t::ANISOTROPIC;
            for(const Json::ArrayIndex i : std::ranges::iota_view{0u, 4u})
                conductivity[i] = conduct[i].template as<T>();
        } else
            throw std::domain_error{"The thermal conductivity should be either a number in the isotropic case, "
                                    "or an array of size 2 in the orthotropic case and size 4 in the anisotropic case"};
    }

    Json::Value save_conductivity() const {
        Json::Value result;
        switch (material) {
        case material_t::ISOTROPIC:
            result = conductivity[0];
        break;

        case material_t::ORTHOTROPIC: {
            result.append(conductivity.front());
            result.append(conductivity.back());
        } break;

        case material_t::ANISOTROPIC:
            for(const T val : conductivity)
                result.append(val);
        break;
        }
        return result;
    }

public:
    material_t material = material_t::ISOTROPIC; // not json field
    std::array<T, 4> conductivity = {T{1}};      // required
    T capacity = T{1};
    T density = T{1};

    explicit constexpr thermal_material_data() noexcept = default;
    explicit thermal_material_data(const Json::Value& physical) {
        check_required_fields(physical, { "conductivity" });
        read_conductivity(physical["conductivity"]);
        capacity = physical.get("capacity", T{1}).template as<T>();
        density = physical.get("density", T{1}).template as<T>();
    }

    Json::Value to_json() const {
        Json::Value result;
        result["conductivity"] = save_conductivity();
        result["capacity"] = capacity;
        result["density"] = density;
        return result;
    }
};

}

#endif
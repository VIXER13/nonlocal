#ifndef NONLOCAL_THERMAL_CONFIG_DATA_HPP
#define NONLOCAL_THERMAL_CONFIG_DATA_HPP

#include "save_data.hpp"
#include "mesh_data.hpp"
#include "time_data.hpp"
#include "boundaries_conditions_data.hpp"
#include "material_data.hpp"
#include "thermal_equation_data.hpp"
#include "thermal_boundary_condition_data.hpp"
#include "thermal_material_data.hpp"

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
using thermal_boundaries_conditions_data = boundaries_conditions_data<thermal_boundary_condition_data, T, Dimension>;

template<std::floating_point T, size_t Dimension>
struct stationary_thermal_data {
    using materials_t = std::conditional_t<
        Dimension == 1,
        std::vector<material_data<thermal_material_data, T, 1>>,
        std::unordered_map<std::string, material_data<thermal_material_data, T, Dimension>>
    >;

    nlohmann::json other;
    save_data save;
    mesh_data<Dimension> mesh;
    thermal_equation_data<T> equation;
    thermal_boundaries_conditions_data<T, Dimension> boundaries; // required
    materials_t materials;                                       // required

    explicit stationary_thermal_data(const nlohmann::json& value)
        : other{value.value("other", nlohmann::json{})}
        , save{value.value("save", nlohmann::json{})}
        , equation{value.value("equation", nlohmann::json{})} {
        if constexpr (Dimension == 1) {
            check_required_fields(value, { "boundaries", "materials" });
            if (value.contains("mesh"))
                mesh = mesh_data<Dimension>{value["mesh"]};
            const nlohmann::json& segments = value["materials"];
            if (!segments.is_array() || segments.empty())
                throw std::domain_error{"Field \"materials\" must be not empty array."};
            materials.reserve(segments.size());
            for(const auto& segment : segments)
                materials.emplace_back(segment);
        } else {
            check_required_fields(value, { "boundaries", "materials", "mesh" });
            mesh = mesh_data<Dimension>{value["mesh"]};
            const nlohmann::json& areas = value["materials"];
            if (!areas.is_object())
                throw std::domain_error{"Field \"materials\" must be a key-value map, where key is material name and value is material parameters."};
            for(const auto& [name, material] : areas.items())
                materials.emplace(name, material);
        }
        boundaries = thermal_boundaries_conditions_data<T, Dimension>{value["boundaries"]};
    }

    virtual ~stationary_thermal_data() noexcept = default;

    operator nlohmann::json() const {
        return {
            {"other", other},
            {"save", save},
            {"mesh", mesh},
            {"equation", equation},
            {"boundaries", boundaries},
            {"materials", materials}
        };
    }
};

template<std::floating_point T, size_t Dimension>
struct nonstationary_thermal_data final : public stationary_thermal_data<T, Dimension> {
    time_data<T> time;

    explicit nonstationary_thermal_data(const nlohmann::json& value)
        : stationary_thermal_data<T, Dimension>{value} {
        check_required_fields(value, { "time" });
        time = time_data<T>{value["time"]};
    }

    ~nonstationary_thermal_data() noexcept override = default;

    operator nlohmann::json() const {
        nlohmann::json result = stationary_thermal_data<T, Dimension>::operator nlohmann::json();
        result["time"] = time;
        return result;
    }
};

}

#endif
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

    Json::Value other;
    save_data save;
    mesh_data<Dimension> mesh;
    thermal_equation_data<T> equation;
    thermal_boundaries_conditions_data<T, Dimension> boundaries; // required
    materials_t materials;                                       // required

    explicit stationary_thermal_data(const Json::Value& value)
        : other{value.get("other", {})}
        , save{value.get("save", {})}
        , equation{value.get("equation", {})} {
        if constexpr (Dimension == 1) {
            check_required_fields(value, { "boundaries", "materials" });
            if (value.isMember("mesh"))
                mesh = mesh_data<Dimension>{value["mesh"]};
            const Json::Value& segments = value["materials"];
            if (!segments.isArray() || segments.empty())
                throw std::domain_error{"Field \"materials\" must be not empty array."};
            materials.reserve(segments.size());
            for(const Json::Value& segment : segments)
                materials.emplace_back(segment);
        } else {
            check_required_fields(value, { "boundaries", "materials", "mesh" });
            mesh = mesh_data<Dimension>{value["mesh"]};
            const Json::Value& areas = value["materials"];
            if (!areas.isObject())
                throw std::domain_error{"Field \"materials\" must be a key-value map, where key is material name and value is material parameters."};
            for(const std::string& name : areas.getMemberNames())
                materials.emplace(name, areas[name]);
        }
        boundaries = thermal_boundaries_conditions_data<T, Dimension>{value["boundaries"]};
    }

    virtual ~stationary_thermal_data() noexcept = default;

    Json::Value to_json() const {
        Json::Value result;
        result["other"] = other;
        result["save"] = save.to_json();
        result["mesh"] = mesh.to_json();
        result["equation"] = equation.to_json();
        result["boundaries"] = boundaries.to_json();
        if constexpr (Dimension == 1) {
            Json::Value& segments = result["materials"] = Json::arrayValue;
            for(const auto& segment : materials)
                segments.append(segment.to_json());
        } else {
            Json::Value& areas = result["materials"] = Json::objectValue;
            for(const auto& [name, area] : materials)
                areas[name] = area.to_json();
        }
        return result;
    }
};

template<std::floating_point T, size_t Dimension>
struct nonstationary_thermal_data final : public stationary_thermal_data<T, Dimension> {
    time_data<T> time;

    explicit nonstationary_thermal_data(const Json::Value& value)
        : stationary_thermal_data<T, Dimension>{value} {
        check_required_fields(value, { "time" });
        time = time_data<T>{value["time"]};
    }

    ~nonstationary_thermal_data() noexcept override = default;

    Json::Value to_json() const {
        Json::Value result = stationary_thermal_data<T, Dimension>::to_json();
        result["time"] = time.to_json();
        return result;
    }
};

}

#endif
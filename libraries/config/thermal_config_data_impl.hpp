#ifndef NONLOCAL_THERMAL_CONFIG_DATA_IMPL_HPP
#define NONLOCAL_THERMAL_CONFIG_DATA_IMPL_HPP

namespace nonlocal::config {

template<std::floating_point T>
thermal_equation_data<T>::thermal_equation_data(const Json::Value& equation) {
    energy = equation.get("energy", T{0}).template as<T>();
    right_part = equation.get("right_part", T{0}).template as<T>();
    initial_distribution = equation.get("initial_distribution", T{0}).template as<T>();
}

template<std::floating_point T>
Json::Value thermal_equation_data<T>::to_json() const {
    Json::Value result;
    result["energy"] = energy;
    result["right_part"] = right_part;
    result["initial_distribution"] = initial_distribution;
    return result;
}

template<std::floating_point T>
thermal_boundary_condition_data<T>::thermal_boundary_condition_data(const Json::Value& condition) {
    check_required_fields(condition, { "kind" });
    switch (kind = get_thermal_condition(condition["kind"])) {
    case thermal::boundary_condition_t::TEMPERATURE: {
        check_required_fields(condition, { "temperature" });
        temperature = condition["temperature"].template as<T>();
    } break;

    case thermal::boundary_condition_t::FLUX: {
        check_required_fields(condition, { "flux" });
        flux = condition["flux"].template as<T>();
    } break;

    case thermal::boundary_condition_t::CONVECTION: {
        check_required_fields(condition, { "temperature", "heat_transfer" });
        temperature = condition["temperature"].template as<T>();
        heat_transfer = condition["heat_transfer"].template as<T>();
    } break;

    case thermal::boundary_condition_t::RADIATION: {
        check_required_fields(condition, { "emissivity" });
        emissivity = condition["emissivity"].template as<T>();
    } break;

    case thermal::boundary_condition_t::COMBINED: {
        temperature = condition.get("temperature", T{0}).template as<T>();
        flux = condition.get("flux", T{0}).template as<T>();
        heat_transfer = condition.get("heat_transfer", T{0}).template as<T>();
        emissivity = condition.get("emissivity", T{0}).template as<T>();
    } break;

    default:
        throw std::domain_error{"Unknown boundary condition type: " + std::to_string(uint(kind))};
    }
}

template<std::floating_point T>
Json::Value thermal_boundary_condition_data<T>::to_json() const {
    Json::Value result;
    result["kind"] = get_thermal_condition(kind);
    result["temperature"] = temperature;
    result["flux"] = flux;
    result["heat_transfer"] = heat_transfer;
    result["emissivity"] = emissivity;
    return result;
}

template<std::floating_point T, size_t Dimension>
thermal_boundaries_conditions_data<T, Dimension>::thermal_boundaries_conditions_data(const Json::Value& boundaries) {
    if constexpr (Dimension == 1)
        check_required_fields(boundaries, { "left", "right" });
    for(const std::string& name : boundaries.getMemberNames())
        conditions[name] = thermal_boundary_condition_data<T>{boundaries[name]};
}

template<std::floating_point T, size_t Dimension>
Json::Value thermal_boundaries_conditions_data<T, Dimension>::to_json() const {
    Json::Value result;
    for(const auto& [name, condition] : conditions)
        result[name] = condition.to_json();
    return result;
}

template<std::floating_point T>
thermal_material_data<T, 1>::thermal_material_data(const Json::Value& physical) {
    check_required_fields(physical, { "conductivity" });
    conductivity = physical["conductivity"].template as<T>();
    capacity = physical.get("capacity", T{1}).template as<T>();
    density = physical.get("density", T{1}).template as<T>();
}

template<std::floating_point T>
Json::Value thermal_material_data<T, 1>::to_json() const {
    Json::Value result;
    result["conductivity"] = conductivity;
    result["capacity"] = capacity;
    result["density"] = density;
    return result;
}

template<std::floating_point T>
void thermal_material_data<T, 2>::read_conductivity(const Json::Value& conduct) {
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

template<std::floating_point T>
Json::Value thermal_material_data<T, 2>::save_conductivity() const {
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

template<std::floating_point T>
thermal_material_data<T, 2>::thermal_material_data(const Json::Value& physical) {
    check_required_fields(physical, { "conductivity" });
    read_conductivity(physical["conductivity"]);
    capacity = physical.get("capacity", T{1}).template as<T>();
    density = physical.get("density", T{1}).template as<T>();
}

template<std::floating_point T>
Json::Value thermal_material_data<T, 2>::to_json() const {
    Json::Value result;
    result["conductivity"] = save_conductivity();
    result["capacity"] = capacity;
    result["density"] = density;
    return result;
}

template<std::floating_point T, size_t Dimension>
stationary_thermal_data<T, Dimension>::stationary_thermal_data(const Json::Value& value)
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

template<std::floating_point T, size_t Dimension>
Json::Value stationary_thermal_data<T, Dimension>::to_json() const {
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

template<std::floating_point T, size_t Dimension>
nonstationary_thermal_data<T, Dimension>::nonstationary_thermal_data(const Json::Value& value)
    : stationary_thermal_data<T, Dimension>{value} {
    check_required_fields(value, { "nonstationary" });
    nonstationary = nonstationary_data<T>{value["nonstationary"]};
}

template<std::floating_point T, size_t Dimension>
Json::Value nonstationary_thermal_data<T, Dimension>::to_json() const {
    Json::Value result = stationary_thermal_data<T, Dimension>::to_json();
    result["nonstationary"] = nonstationary.to_json();
    return result;
}

}

#endif
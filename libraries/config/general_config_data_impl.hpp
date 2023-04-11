#ifndef NONLOCAL_GENERAL_CONFIG_DATA_IMPL_HPP
#define NONLOCAL_GENERAL_CONFIG_DATA_IMPL_HPP

namespace nonlocal::config {

template<size_t Dimension>
mesh_data<Dimension>::mesh_data(const Json::Value& value) {
    check_required_fields(value, { "mesh" });
    mesh = value["mesh"].asString();
}

template<size_t Dimension>
Json::Value mesh_data<Dimension>::to_json() const {
    Json::Value result;
    result["mesh"] = mesh.string();
    return result;
}

template<>
struct mesh_data<1u> final {
    size_t element_order = 1;
    size_t quadrature_order = 1;

    explicit constexpr mesh_data() noexcept = default;
    explicit mesh_data(const Json::Value& value)
        : element_order{get_order(value.get("element_order", 1))}
        , quadrature_order{get_order(value.get("quadrature_order", element_order))} {}

    Json::Value to_json() const {
        Json::Value result;
        result["element_order"] = get_order(element_order);
        result["quadrature_order"] = get_order(quadrature_order);
        return result;
    }
};

template<template<class, size_t> class Condition, std::floating_point T, size_t Dimension>
boundaries_conditions_data<Condition, T, Dimension>::boundaries_conditions_data(const Json::Value& boundaries) {
    if constexpr (Dimension == 1)
        check_required_fields(boundaries, { "left", "right" });
    for(const std::string& name : boundaries.getMemberNames())
        conditions[name] = Condition<T, Dimension>{boundaries[name]};
}

template<template<class, size_t> class Condition, std::floating_point T, size_t Dimension>
Json::Value boundaries_conditions_data<Condition, T, Dimension>::to_json() const {
    Json::Value result;
    for(const auto& [name, condition] : conditions)
        result[name] = condition.to_json();
    return result;
}

template<std::floating_point T>
time_data<T>::time_data(const Json::Value& nonstationary) {
    check_required_fields(nonstationary, { "time_step", "steps_count"});
    time_step = nonstationary["time_step"].template as<T>();
    initial_time = nonstationary.get("initial_time", T{0}).template as<T>();
    steps_count = nonstationary["steps_count"].asUInt64();
    save_frequency = nonstationary.get("save_frequency", 1u).asUInt64();
}

template<std::floating_point T>
Json::Value time_data<T>::to_json() const {
    Json::Value result;
    result["time_step"] = time_step;
    result["initial_time"] = initial_time;
    result["steps_count"] = steps_count;
    result["save_frequency"] = save_frequency;
    return result;
}

template<std::floating_point T, size_t Dimension>
model_data<T, Dimension>::radius_t model_data<T, Dimension>::read_radius(const Json::Value& arr, const std::string& field) {
    std::array<T, Dimension> result;
    
    if (arr.isDouble())
        result.fill(arr.template as<T>());
    else if (arr.isArray() && arr.size() == Dimension)
        for(const Json::ArrayIndex i : std::ranges::iota_view{0u, Dimension})
            result[i] = arr[i].template as<T>();
    else
        throw std::domain_error{"Field \"" + field + "\" must be an array with length " + std::to_string(Dimension)};

    if constexpr (std::is_same_v<radius_t, T>)
        return result.front();
    else
        return result;
}

template<std::floating_point T, size_t Dimension>
model_data<T, Dimension>::model_data(const Json::Value& model) {
    check_required_fields(model, { "local_weight", "nonlocal_radius" });
    local_weight = model["local_weight"].template as<T>();
    nonlocal_radius = read_radius(model["nonlocal_radius"], "nonlocal_radius");
    search_radius = !model.isMember("search_radius") ? nonlocal_radius :
                    read_radius(model["search_radius"], "search_radius");
}

template<std::floating_point T, size_t Dimension>
Json::Value model_data<T, Dimension>::to_json() const {
    Json::Value result;
    result["local_weight"] = local_weight;
    if constexpr (std::is_same_v<radius_t, T>) {
        result["nonlocal_radius"] = nonlocal_radius;
        result["search_radius"] = search_radius;
    } else {
        Json::Value& nonloc = result["nonlocal_radius"] = Json::arrayValue;
        Json::Value& search = result["search_radius"] = Json::arrayValue;
        for(const size_t i : std::ranges::iota_view{0u, Dimension}) {
            nonloc.append(nonlocal_radius[i]);
            search.append(search_radius[i]);
        }
    }
    return result;
}

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
material_data<Physics, T, Dimension>::material_data(const Json::Value& material) {
    check_required_fields(material, { "physical" });
    physical = Physics<T, Dimension>{material["physical"]};
    if (material.isMember("model"))
        model = model_data<T, Dimension>{material["model"]};
}

template<template<class, size_t> class Physics, std::floating_point T, size_t Dimension>
Json::Value material_data<Physics, T, Dimension>::to_json() const {
    Json::Value result;
    result["physical"] = physical.to_json();
    result["model"] = model.to_json();
    return result;
}

template<template<class, size_t> class Physics, std::floating_point T>
material_data<Physics, T, 1>::material_data(const Json::Value& material) {
    check_required_fields(material, { "elements_count", "length", "physical" });
    elements_count = material["elements_count"].asUInt64();
    length = material["length"].template as<T>();
    physical = Physics<T, 1>{material["physical"]};
    if (material.isMember("model"))
        model = model_data<T, 1>{material["model"]};
}

template<template<class, size_t> class Physics, std::floating_point T>
Json::Value material_data<Physics, T, 1>::to_json() const {
    Json::Value result;
    result["elements_count"] = elements_count;
    result["length"] = length;
    result["physical"] = physical.to_json();
    result["model"] = model.to_json();
    return result;
}

}

#endif
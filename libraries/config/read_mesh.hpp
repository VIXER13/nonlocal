#pragma once

#include "config_utils.hpp"

#include <metamath/metamath.hpp>
#include <logger/logger.hpp>

#include <mesh/mesh_1d/mesh_1d.hpp>
#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::config {

enum class order_t : uint8_t {
    Unknown,
    Linear,
    Quadratic,
    Cubic,
    Quartic,
    Quintic
};

NLOHMANN_JSON_SERIALIZE_ENUM(order_t, {
    {order_t::Unknown, nullptr},
    {order_t::Linear, "linear"},
    {order_t::Quadratic, "quadratic"},
    {order_t::Cubic, "cubic"},
    {order_t::Quartic, "quartic"},
    {order_t::Quintic, "quintic"}
})

class _read_mesh_1d final {
    template<class T>
    using finite_element_1d = metamath::finite_element::element_1d_integrate<T>;

    template<class T, size_t Order>
    using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, Order>;
    template<class T>
    using element_1d_integrate = metamath::finite_element::element_1d_integrate<T>;
    template<class T, size_t Order>
    using element_1d = metamath::finite_element::element_1d<T, metamath::finite_element::lagrangian_element_1d, Order>;

    static bool is_valid_order(const size_t order) noexcept;
    static order_t get_order(const nlohmann::json& config, const std::string& field, const order_t default_order = order_t::Linear);

    template<std::floating_point T>
    static finite_element_1d<T> read_element(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    static T read_search_radius(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static mesh::segment_data<T> read_segment(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    static std::vector<mesh::segment_data<T>> read_segments(const nlohmann::json& config, const std::string& path);

    explicit constexpr _read_mesh_1d() noexcept = default;
    
public:
    template<std::floating_point T>
    friend std::shared_ptr<mesh::mesh_1d<T>> read_mesh_1d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
_read_mesh_1d::finite_element_1d<T> _read_mesh_1d::read_element(const nlohmann::json& config, const std::string& path) {
    if (!config.empty())
        check_optional_fields(config, {"element_order", "quadrature_order"}, append_access_sign(path));
    const size_t element_order = size_t(_read_mesh_1d::get_order(config, "element_order"));
    const size_t quadrature_order = size_t(_read_mesh_1d::get_order(config, "quadrature_order", order_t(element_order)));
    if (quadrature_order < element_order)
        logger::warning() << "The order of the quadrature is lower than the order of the element, this may lead to computational problems." << std::endl;
    return metamath::finite_element::make_element_1d_integrated<T>(element_order, quadrature_order);
}

template<std::floating_point T>
T _read_mesh_1d::read_search_radius(const nlohmann::json& config, const std::string& path) {
    if (config.contains("search_radius"))
        return config["search_radius"].get<T>();
    if (config.contains("nonlocal_radius"))
        return config["nonlocal_radius"].get<T>();
    return 0;
}

template<std::floating_point T>
mesh::segment_data<T> _read_mesh_1d::read_segment(const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, {"elements_count", "length"}, path);
    check_optional_fields(config, {"model"}, path);
    return {
        .length = config["length"].get<T>(),
        .search_radius = read_search_radius<T>(config.value("model", nlohmann::json::object()), "model"),
        .elements = config["elements_count"].get<size_t>()
    };
}

template<std::floating_point T>
std::vector<mesh::segment_data<T>> _read_mesh_1d::read_segments(const nlohmann::json& config, const std::string& path) {
    if (!config.is_array() || config.empty())
        throw std::domain_error{"segments initialization requires the initializing config to be an nonempty array."};
    std::vector<mesh::segment_data<T>> segments;
    segments.reserve(config.size());
    for(const size_t i : std::ranges::iota_view{0u, segments.capacity()})
        segments.push_back(read_segment<T>(config[i], append_access_sign(path, i)));
    return segments;
}

template<std::floating_point T>
std::shared_ptr<mesh::mesh_1d<T>> read_mesh_1d(const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, {"materials"});
    check_optional_fields(config, {"mesh"});
    const std::string path_with_access = append_access_sign(path);
    return std::make_shared<mesh::mesh_1d<T>>(
        _read_mesh_1d::read_element<T>(config.value("mesh", nlohmann::json::object()), path_with_access + "mesh"),
        _read_mesh_1d::read_segments<T>(config["materials"], path_with_access + "materials")
    );
}

template<std::floating_point T, std::integral I>
std::shared_ptr<mesh::mesh_2d<T, I>> read_mesh_2d(const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, { "path" }, append_access_sign(path));
    return std::make_shared<mesh::mesh_2d<T, I>>(config["path"].get<std::string>());
}

}
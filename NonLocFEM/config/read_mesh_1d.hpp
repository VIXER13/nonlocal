#pragma once

#include "config_utils.hpp"

#include "logger.hpp"
#include "mesh_1d.hpp"
#include "metamath.hpp"

namespace nonlocal::config {

enum class order_t : uint8_t {
    Unknown,
    Linear,
    Quadratic,
    Сubic,
    Quartic,
    Quintic
};

NLOHMANN_JSON_SERIALIZE_ENUM(order_t, {
    {order_t::Unknown, nullptr},
    {order_t::Linear, "linear"},
    {order_t::Quadratic, "quadratic"},
    {order_t::Сubic, "сubic"},
    {order_t::Quartic, "quartic"},
    {order_t::Quintic, "quintic"}
})

template<class T>
using quadrature_1d_ptr = std::unique_ptr<metamath::finite_element::quadrature_1d_base<T>>;
template<class T>
using finite_element_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

class _read_mesh_1d final {
    template<class T, size_t Order>
    using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, Order>;
    template<class T, size_t Order>
    using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, Order>;

    explicit constexpr _read_mesh_1d() noexcept = default;

    static bool is_valid_order(const size_t order) noexcept;
    static order_t get_order(const nlohmann::json& config, const std::string& field, const order_t default_order = order_t::Linear);

    template<std::floating_point T>
    friend quadrature_1d_ptr<T> make_quadrature(const order_t order);
    template<std::floating_point T>
    friend finite_element_1d_ptr<T> make_element(const order_t order, const quadrature_1d_ptr<T>& quadrature);
    template<std::floating_point T>
    friend finite_element_1d_ptr<T> make_element(const order_t element_order, const order_t quadrature_order);

    template<std::floating_point T>
    friend finite_element_1d_ptr<T> read_element(const nlohmann::json& config, const std::string& path);

    template<std::floating_point T>
    friend T read_search_radius(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    friend mesh::segment_data<T> read_segment(const nlohmann::json& config, const std::string& path);
    template<std::floating_point T>
    friend std::vector<mesh::segment_data<T>> read_segments(const nlohmann::json& config, const std::string& path);
    
public:
    template<std::floating_point T>
    friend std::shared_ptr<mesh::mesh_1d<T>> read_mesh_1d(const nlohmann::json& config, const std::string& path);
};

template<std::floating_point T>
quadrature_1d_ptr<T> make_quadrature(const order_t order) {
    switch(order) {
    case order_t::Linear:
        return std::make_unique<_read_mesh_1d::quadrature<T, 1>>();
    case order_t::Quadratic:
        return std::make_unique<_read_mesh_1d::quadrature<T, 2>>();
    case order_t::Сubic:
        return std::make_unique<_read_mesh_1d::quadrature<T, 3>>();
    case order_t::Quartic:
        return std::make_unique<_read_mesh_1d::quadrature<T, 4>>();
    case order_t::Quintic:
        return std::make_unique<_read_mesh_1d::quadrature<T, 5>>();
    }
    throw std::logic_error{"Invalid quadrature order: " + std::to_string(size_t(order))};
}

template<std::floating_point T>
finite_element_1d_ptr<T> make_element(const order_t order, const quadrature_1d_ptr<T>& quadrature) {
    switch(order) {
    case order_t::Linear:
        return std::make_unique<_read_mesh_1d::element_1d<T, 1>>(*quadrature);
    case order_t::Quadratic:
        return std::make_unique<_read_mesh_1d::element_1d<T, 2>>(*quadrature);
    case order_t::Сubic:
        return std::make_unique<_read_mesh_1d::element_1d<T, 3>>(*quadrature);
    case order_t::Quartic:
        return std::make_unique<_read_mesh_1d::element_1d<T, 4>>(*quadrature);
    case order_t::Quintic:
        return std::make_unique<_read_mesh_1d::element_1d<T, 5>>(*quadrature);
    }
    throw std::logic_error{"Invalid element order: " + std::to_string(size_t(order))};
}

template<std::floating_point T>
finite_element_1d_ptr<T> make_element(const order_t element_order, const order_t quadrature_order) {
    return make_element(element_order, make_quadrature<T>(quadrature_order == order_t::Unknown ? element_order : quadrature_order));
}

template<std::floating_point T>
finite_element_1d_ptr<T> read_element(const nlohmann::json& config, const std::string& path) {
    if (!config.empty())
        check_optional_fields(config, {"element_order", "quadrature_order"}, append_access_sign(path));
    const order_t element_order = _read_mesh_1d::get_order(config, "element_order");
    const order_t quadrature_order = _read_mesh_1d::get_order(config, "quadrature_order", element_order);
    if (size_t(quadrature_order) < size_t(element_order))
        logger::warning() << "The order of the quadrature is lower than the order of the element, this may lead to computational problems." << std::endl;
    return make_element<T>(element_order, quadrature_order);
}

template<std::floating_point T>
T read_search_radius(const nlohmann::json& config, const std::string& path) {
    if (config.contains("search_radius"))
        return config["search_radius"].get<T>();
    if (config.contains("nonlocal_radius"))
        return config["nonlocal_radius"].get<T>();
    return 0;
}

template<std::floating_point T>
mesh::segment_data<T> read_segment(const nlohmann::json& config, const std::string& path) {
    check_required_fields(config, {"elements_count", "length"}, path);
    check_optional_fields(config, {"model"}, path);
    return {
        .length = config["length"].get<T>(),
        .search_radius = read_search_radius<T>(config.value("model", nlohmann::json::object()), "model"),
        .elements = config["elements_count"].get<size_t>()
    };
}

template<std::floating_point T>
std::vector<mesh::segment_data<T>> read_segments(const nlohmann::json& config, const std::string& path) {
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
        read_element<T>(config.value("mesh", nlohmann::json::object()), path_with_access + "mesh"),
        read_segments<T>(config["materials"], path_with_access + "materials")
    );
}

}
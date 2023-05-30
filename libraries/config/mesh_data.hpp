#ifndef NONLOCAL_CONFIG_MESH_DATA_HPP
#define NONLOCAL_CONFIG_MESH_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

template<size_t Dimension>
struct mesh_data final {
    std::filesystem::path path; // required

    explicit mesh_data() = default;
    explicit mesh_data(const nlohmann::json& value) {
        check_required_fields(value, { "path" });
        path = value["path"].get<std::string>();
    }

    operator nlohmann::json() const {
        return { {"path", path.string()} };
    }
};

template<>
struct mesh_data<1u> final {
    size_t element_order = 1;
    size_t quadrature_order = 1;

    explicit constexpr mesh_data() noexcept = default;
    explicit mesh_data(const nlohmann::json& value)
        : element_order{get_order(value.value<nlohmann::json>("element_order", 1))}
        , quadrature_order{get_order(value.value<nlohmann::json>("quadrature_order", element_order))} {}

    operator nlohmann::json() const {
        return {
            {"element_order", get_order(element_order)},
            {"quadrature_order", get_order(quadrature_order)}
        };
    }
};

}

#endif
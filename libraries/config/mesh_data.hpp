#ifndef NONLOCAL_CONFIG_MESH_DATA_HPP
#define NONLOCAL_CONFIG_MESH_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

enum class order_t : uint8_t {
    UNKNOWN,
    LINEAR,
    QUADRATIC,
    QUBIC,
    QUARTIC,
    QUINTIC
};

NLOHMANN_JSON_SERIALIZE_ENUM(order_t, {
    {order_t::UNKNOWN, nullptr},
    {order_t::LINEAR, "linear"},
    {order_t::QUADRATIC, "quadratic"},
    {order_t::QUBIC, "qubic"},
    {order_t::QUARTIC, "quartic"},
    {order_t::QUINTIC, "quintic"}
})

template<size_t Dimension>
struct mesh_data final {
    std::filesystem::path path; // required

    explicit mesh_data() = default;
    explicit mesh_data(const nlohmann::json& config, const std::string& path = "") {
        check_required_fields(config, { "path" }, append_access_sign(path));
        path = config["path"].get<std::string>();
    }

    operator nlohmann::json() const {
        return { {"path", path.string()} };
    }
};

template<>
struct mesh_data<1u> final {
    order_t element_order = order_t::LINEAR;
    order_t quadrature_order = element_order;

    static bool is_valid_order(const size_t order) noexcept;

    explicit constexpr mesh_data() noexcept = default;
    explicit mesh_data(const nlohmann::json& config, const std::string& path = "");

    operator nlohmann::json() const;
};

}

#endif
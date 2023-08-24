#include "mesh_data.hpp"

namespace {

using namespace nonlocal::config;

order_t get_order(const nlohmann::json& config, const std::string& field, const order_t default_order = order_t::LINEAR) {
    if (config.contains("element_order"))
        if (const nlohmann::json& value = config["element_order"]; value.is_number_integer())
            if (const size_t order = value.get<size_t>(); mesh_data<1u>::is_valid_order(order))
                return order_t(order);
    return default_order;
}

}

namespace nonlocal::config {

bool mesh_data<1u>::is_valid_order(const size_t order) noexcept {
    return order > 0 && order < 6;
}

mesh_data<1u>::mesh_data(const nlohmann::json& config, const std::string& path)
    : element_order{get_order(config, "element_order")}
    , quadrature_order{get_order(config, "quadrature_order", element_order)} {
    check_optional_fields(config, {"element_order", "quadrature_order"}, append_access_sign(path));
}

mesh_data<1u>::operator nlohmann::json() const {
    return {
        {"element_order", element_order},
        {"quadrature_order", quadrature_order}
    };
}

}
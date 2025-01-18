#include "read_mesh_1d.hpp"

namespace nonlocal::config {

bool _read_mesh_1d::is_valid_order(const size_t order) noexcept {
    return order > 0 && order < 6;
}

order_t _read_mesh_1d::get_order(const nlohmann::json& config, const std::string& field, const order_t default_order) {
    if (config.contains(field)) {
        if (const nlohmann::json& value = config[field]; value.is_number_integer()) {
            if (const size_t order = value.get<size_t>(); is_valid_order(order))
                return order_t(order);
            else
                throw std::domain_error{"Unsupported \"" + field + "\": " + std::to_string(order)};
        } else if (value.is_string()) {
            if(const order_t order = value.get<order_t>(); order != order_t::Unknown)
                return order;
            else
                throw std::domain_error{"Unsupported \"" + field + "\": " + value.get<std::string>()};
        }
        throw std::domain_error{"Unsupported \"" + field + "\" input."};
    }
    return default_order;
}

}
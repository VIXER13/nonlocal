#include "mesh_data.hpp"

namespace nonlocal::config {

mesh_data<1u>::mesh_data(const nlohmann::json& value)
    : element_order{value.value<nlohmann::json>("element_order", order_t::LINEAR)}
    , quadrature_order{value.value<nlohmann::json>("quadrature_order", element_order)} {}

mesh_data<1u>::operator nlohmann::json() const {
    return {
        {"element_order", element_order},
        {"quadrature_order", quadrature_order}
    };
}

}
#ifndef MECHANICAL_BOUNDARY_CONDITION_DATA_HPP
#define MECHANICAL_BOUNDARY_CONDITION_DATA_HPP

#include "config_utils.hpp"

namespace nonlocal::config {

enum class mechanical_boundary_condition_t : uint8_t {
    UNKNOWN,
    DISPLACEMENT,
    PRESSURE
};

NLOHMANN_JSON_SERIALIZE_ENUM(mechanical_boundary_condition_t, {
    {mechanical_boundary_condition_t::UNKNOWN, nullptr},
    {mechanical_boundary_condition_t::DISPLACEMENT, "displacement"},
    {mechanical_boundary_condition_t::PRESSURE, "pressure"}
})

}

#endif
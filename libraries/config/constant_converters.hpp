#ifndef NONLOCAL_CONFIG_CONSTANT_CONVERTERS_HPP
#define NONLOCAL_CONFIG_CONSTANT_CONVERTERS_HPP

#include "nonlocal_constants.hpp"

#include <nlohmann/json.hpp>

namespace nonlocal::thermal {

NLOHMANN_JSON_SERIALIZE_ENUM(thermal::boundary_condition_t, {
    {TEMPERATURE, "temperature"},
    {FLUX, "flux"},
    {CONVECTION, "convection"},
    {RADIATION, "radiation"},
    {COMBINED, "combined"}
})

}

#endif
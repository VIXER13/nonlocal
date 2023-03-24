#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

#include "metamath.hpp"
#include "nonlocal_constants.hpp"
#include "../../equation_parameters.hpp"

#include <string>
#include <unordered_map>

namespace nonlocal::thermal {

template<class T>
struct parameter_2d final {
    metamath::types::square_matrix<T, 2> conductivity = {T{1}};
    T capacity = T{1};
    T density = T{1};
    material_t material = material_t::ISOTROPIC;
};

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

};

#endif
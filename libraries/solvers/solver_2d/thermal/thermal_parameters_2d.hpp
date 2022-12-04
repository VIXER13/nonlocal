#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

#include "metamath.hpp"
#include "../../solvers_constants.hpp"


namespace nonlocal::thermal {

template<class T>
struct parameter_2d final {
    metamath::types::square_matrix<T, 2> conductivity = {T{1}};
    T capacity = T{1};
    T density = T{1};
    material_t material = material_t::ISOTROPIC;
};

};

#endif
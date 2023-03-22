#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

#include "../../equation_parameters.hpp"

namespace nonlocal::thermal {

template<class T>
struct parameter_1d final {
    T conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
};

template<class T>
using parameters_1d = std::vector<equation_parameters<1, T, parameter_1d>>;

}

#endif
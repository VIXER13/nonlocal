#ifndef NONLOCAL_THERMAL_PARAMETERS_1D
#define NONLOCAL_THERMAL_PARAMETERS_1D

namespace nonlocal::thermal {

template<class T>
struct parameters_1d final {
    T conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
};


}

#endif
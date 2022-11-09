#ifndef HEAT_EQUATION_PARAMETERS_1D_HPP
#define HEAT_EQUATION_PARAMETERS_1D_HPP

#include <array>

namespace nonlocal::thermal {

template<class T>
struct nonstationary_equation_parameters_1d final {
    T conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
    std::array<T, 2> transfer = {T{0}, T{0}};
};

}

#endif
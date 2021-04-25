#ifndef HEAT_EQUATION_PARAMETERS_HPP
#define HEAT_EQUATION_PARAMETERS_HPP

#include <array>

namespace nonlocal::heat {

enum class calc_type : bool { ISOTROPIC, ORTHOTROPIC };

template<class T, calc_type Calc_Type>
struct calculation_parameters final {
    std::array<T, Calc_Type == calc_type::ORTHOTROPIC ? 2 : 1> lambda;
    T rho = 1,
      c   = 1;

    calculation_parameters() {
        lambda.fill(T{1});
    }
};

}

#endif
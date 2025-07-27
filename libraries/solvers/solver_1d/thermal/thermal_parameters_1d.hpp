#pragma once

#include <solvers/base/equation_parameters.hpp>

#include <memory>

namespace nonlocal::thermal {

template<std::floating_point T>
struct parameter_1d final {
    static constexpr size_t Dimension = 1;
    coefficient_t<T, Dimension> conductivity;
    T capacity = T{1};
    T density = T{1};
    T relaxation_time = T{0};
};

template<class T>
using parameters_1d = std::vector<equation_parameters<1, T, parameter_1d>>;

}
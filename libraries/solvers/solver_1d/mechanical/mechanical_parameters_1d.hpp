#pragma once

#include <solvers/base/equation_parameters.hpp>

#include <memory>

namespace nonlocal::solver_1d::mechanical {

template<std::floating_point T>
struct parameter_1d final {
    static constexpr size_t Dimension = 1;
    coefficient_t<T, Dimension> youngs_modulus;
    coefficient_t<T, Dimension> density;
};

template<class T>
using parameters_1d = std::vector<equation_parameters<1, T, parameter_1d>>;

}
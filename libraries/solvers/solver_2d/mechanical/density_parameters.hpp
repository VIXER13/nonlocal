#pragma once

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::mechanical {

template<std::floating_point T>
using raw_density_t = std::variant<
    std::monostate,
    coefficient_t<T, 2>
>;

template<std::floating_point T>
using evaluated_density_t = std::variant<
    std::monostate,
    evaluated_parameters<T>
>;

}
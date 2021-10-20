#ifndef BOUNDARY_CONDITION_1D_HPP
#define BOUNDARY_CONDITION_1D_HPP

#include "../solvers_constants.hpp"
#include <utility>
#include <functional>

namespace nonlocal {

using stationary_boundary_1d = std::pair<boundary_condition_t, T>;
using nonstatinary_boundary_t = std::pair<boundary_condition_t, std::function<T(T)>;

}

#endif
#pragma once

#include <metamath/metamath.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>

namespace unit_tests {

// Mesh initialization is separated into a separate translation unit to avoid the SIOF problem.
template<std::floating_point T, std::signed_integral I>
std::shared_ptr<nonlocal::mesh::mesh_2d<T, I>> init_2d_mesh(std::stringstream& stream, const nonlocal::mesh::mesh_format format);

}
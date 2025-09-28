#pragma once

#include <logger/logger.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>

namespace nonlocal {
    
template<std::floating_point T, std::signed_integral I>
solver_2d::mechanical::mechanical_solution_2d<T> solve_mechanical_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const solver_2d::mechanical::mechanical_parameters_2d<T>& parameters,
    const solver_2d::mechanical::mechanical_boundaries_conditions_2d<T>& boundaries_conditions) {
    auto solution = solver_2d::mechanical::equilibrium_equation<I>(mesh, parameters, boundaries_conditions);
    solution.calc_strain_and_stress();
    return solution;
}
    
}
#pragma once

#include <logger/logger.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>

namespace nonlocal::mechanical {
    
template<std::floating_point T, std::signed_integral I>
mechanical::mechanical_solution_2d<T> solve_mechanical_2d_problem(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh,
    const mechanical_parameters_2d<T>& parameters,
    const mechanical::mechanical_boundaries_conditions_2d<T>& boundaries_conditions) {
    auto solution = nonlocal::mechanical::equilibrium_equation<I>(
        mesh, parameters, boundaries_conditions,
        [](const std::array<T, 2>&) constexpr noexcept { return std::array<T, 2>{}; }
    );
    solution.calc_strain_and_stress();
    return solution;
}
    
}
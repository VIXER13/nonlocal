#ifndef NONLOCAL_STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define NONLOCAL_STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "solvers_utils.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "thermal_boundary_conditions_2d.hpp"
#include "boundary_condition_first_kind_2d.hpp"
#include "boundary_condition_second_kind_2d.hpp"
#include "convection_condition_2d.hpp"
#include "right_part_2d.hpp"
#include "heat_equation_solution_2d.hpp"

#include "conjugate_gradient.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class Matrix_Index, class T, class I, class Right_Part>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                   const parameters_2d<T>& parameters,
                                                                   const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                   const Right_Part& right_part,
                                                                   const T energy = T{0}) {
    static constexpr size_t DoF = 1;
    static constexpr auto is_second_kind = [](const auto& condition) {
        return bool(dynamic_cast<const flux_2d<T>*>(condition.get()));
    };
    const auto conditions = boundaries_conditions | std::views::values;
    const bool is_neumann = std::all_of(conditions.begin(), conditions.end(), is_second_kind);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count() + is_neumann);
    boundary_condition_second_kind_2d(f, *mesh, boundaries_conditions);
    if (is_neumann) {
    //     if (!is_solvable_neumann_problem(*mesh_proxy, f))
    //         throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
        f[f.size() - 1] = energy;
    }

    static constexpr auto check_nonlinear = [](const auto& parameter) { return parameter.second.physical->type != coefficients_t::CONSTANTS; };
    const bool is_nonlinear = std::any_of(parameters.begin(), parameters.end(), check_nonlinear);
    static constexpr auto check_nonlocal = [](const auto& theory) noexcept { return theory.second == theory_t::NONLOCAL; };
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    const bool is_nonlocal = std::any_of(theories.begin(), theories.end(), check_nonlocal);
    const bool is_symmetric = !(is_nonlinear && is_nonlocal);

    auto start_time = std::chrono::high_resolution_clock::now();
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh};
    conductivity.compute(parameters, utils::inner_nodes(mesh->container(), boundaries_conditions), is_symmetric, is_neumann);
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Conductivity matrix calculated time: " << elapsed_seconds.count() << 's' << std::endl;
    convection_condition_2d(conductivity.matrix_inner(), *mesh, boundaries_conditions);
    integrate_right_part<DoF>(f, *mesh, right_part);
    if (!is_neumann)
        boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, conductivity.matrix_bound());

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature;
    start_time = std::chrono::high_resolution_clock::now();
    if (is_symmetric) {
        std::cout << "symmetric problem" << std::endl;
        const slae::conjugate_gradient<T, Matrix_Index> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
        std::cout << "Iterations: " << solver.iterations() << std::endl;
    } else {
        std::cout << "asymmetric problem" << std::endl;
        const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, I>> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
        std::cout << "Iterations: " << solver.iterations() << std::endl;
    }
    elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "SLAE solution time: " << elapsed_seconds.count() << 's' << std::endl;
    return heat_equation_solution_2d<T, I>{mesh, parameters, temperature};
}

}

#endif
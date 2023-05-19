#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "mesh_1d_utils.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_1d final {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T energy = T{0};
};

template<class T, class I>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const parameters_1d<T>& parameters,
                                                                const thermal_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                                const stationary_equation_parameters_1d<T>& additional_parameters) {
    static constexpr auto is_flux = [](const auto& condition) noexcept {
        return bool(dynamic_cast<const flux_1d<T>*>(condition.get()));
    };
    const bool is_neumann = std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_flux);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(f, boundaries_conditions, is_neumann);
    if (additional_parameters.right_part)
        integrate_right_part(f, *mesh, *additional_parameters.right_part);
    if (is_neumann && std::abs(std::reduce(f.begin(), f.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};
    if (is_neumann) 
        f[f.size() - 1] = additional_parameters.energy;
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_f = f;

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(f.size());
    if (additional_parameters.initial_distribution)
        for(const size_t node : mesh->nodes())
            temperature_prev[node] = (*additional_parameters.initial_distribution)(mesh->node_coord(node));
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = temperature_prev;

    static constexpr auto is_solution_depend_parameter = [](const auto& parameter) noexcept {
        return parameter.physical->type == coefficients_t::SOLUTION_DEPENDENT;
    };
    const bool is_sol_depend = std::any_of(parameters.begin(), parameters.end(), is_solution_depend_parameter);

    static constexpr auto check_nonlinear = [](const auto& parameter) { return parameter.physical->type != coefficients_t::CONSTANTS; };
    const bool is_nonlinear = std::any_of(parameters.begin(), parameters.end(), check_nonlinear);
    static constexpr auto check_nonlocal = [](const theory_t theory) noexcept { return theory == theory_t::NONLOCAL; };
    const std::vector<theory_t> theories = theories_types(parameters);
    const bool is_nonlocal = std::any_of(theories.begin(), theories.end(), check_nonlocal);
    const bool is_symmetric = !(is_nonlinear && is_nonlocal);

    T difference = T{1};
    size_t iteration = 0;
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    auto start_time = std::chrono::high_resolution_clock::now();
    while (iteration < additional_parameters.max_iterations && 
           difference > additional_parameters.tolerance) {
        std::swap(temperature_prev, temperature_curr);
        std::copy(initial_f.begin(), initial_f.end(), f.begin());
        using namespace nonlocal::mesh::utils;
        conductivity.template calc_matrix(
            parameters,
            { bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.front().get())),
              bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.back ().get())) },
            is_neumann, is_symmetric,
            is_sol_depend ? std::optional{from_nodes_to_qnodes(*mesh, temperature_prev)} : std::nullopt
        );

        if (is_neumann) {
            std::cout << "neumann problem" << std::endl;
            if (is_symmetric) {
                std::cout << "symmetric problem" << std::endl;
                const Eigen::ConjugateGradient<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>,
                    Eigen::Upper
                > solver{conductivity.matrix_inner()};
                temperature_curr = solver.solveWithGuess(f, temperature_prev);
            } else {
                std::cout << "asymmetric problem" << std::endl;
                const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, I>> solver{conductivity.matrix_inner()};
                temperature_curr = solver.solveWithGuess(f, temperature_prev);
            }
        } else {
            convection_condition_1d(conductivity.matrix_inner(), boundaries_conditions);
            boundary_condition_first_kind_1d(f, conductivity.matrix_bound(), boundaries_conditions);
            if (is_symmetric) {
                const Eigen::SimplicialCholesky<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>, 
                    Eigen::Upper, 
                    Eigen::NaturalOrdering<I>
                > solver{conductivity.matrix_inner()};
                temperature_curr = solver.solve(f);
            } else {
                const Eigen::SparseLU<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>,
                    Eigen::NaturalOrdering<I>
                > solver{conductivity.matrix_inner()};
                temperature_curr = solver.solve(f);
            }
        }

        if (!is_sol_depend)
            break;

        ++iteration;
        difference = (temperature_curr - temperature_prev).norm() / (temperature_curr.norm() ?: T{1});
        std::cout << "norm(prev - curr) = " << difference << std::endl;
    }
    std::cout << "Iterations: " << iteration << std::endl;
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Time: " << elapsed_seconds.count() << 's' << std::endl;
    return heat_equation_solution_1d<T>{mesh, parameters, temperature_curr};
}

}

#endif
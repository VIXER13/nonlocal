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

template<class T>
struct stationary_equation_parameters_2d final {
    std::optional<std::function<T(const std::array<T, 2>&)>> right_part;
    std::optional<std::function<T(const std::array<T, 2>&)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T energy = T{0};
};

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


template<class Matrix_Index, class T, class I>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_nonlinear_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                             const parameters_2d<T>& parameters,
                                                                             const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                                                                             const stationary_equation_parameters_2d<T>& additional_parameters) {
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
        f[f.size() - 1] = additional_parameters.energy;
    }
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_f = f;

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(f.size());
    if (additional_parameters.initial_distribution)
        for(const size_t node : mesh->container().nodes()) 
            temperature_prev[node] = (*(additional_parameters.initial_distribution))(mesh->container().node_coord(node));
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = temperature_prev;

    static constexpr auto is_solution_depend_parameter = [](const auto& parameter) noexcept {
        return parameter.second.physical->type == coefficients_t::SOLUTION_DEPENDENT;
    };
    const bool is_sol_depend = std::any_of(parameters.begin(), parameters.end(), is_solution_depend_parameter);

    static constexpr auto check_nonlinear = [](const auto& parameter) { return parameter.second.physical->type != coefficients_t::CONSTANTS; };
    const bool is_nonlinear = std::any_of(parameters.begin(), parameters.end(), check_nonlinear);
    static constexpr auto check_nonlocal = [](const auto& theory) noexcept { return theory.second == theory_t::NONLOCAL; };
    const std::unordered_map<std::string, theory_t> theories = theories_types(parameters);
    const bool is_nonlocal = std::any_of(theories.begin(), theories.end(), check_nonlocal);
    const bool is_symmetric = !(is_nonlinear && is_nonlocal);

    integrate_right_part<DoF>(f, *mesh, *(additional_parameters.right_part));


    T difference = T{1};
    size_t iteration = 0;
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh};

    auto start_time = std::chrono::high_resolution_clock::now();
    while (iteration < additional_parameters.max_iterations && 
           difference > additional_parameters.tolerance) {
        std::cout << " --------------" << "Iteration №" << iteration << " --------------" << std::endl;
        std::swap(temperature_prev, temperature_curr);
        std::copy(initial_f.begin(), initial_f.end(), f.begin());
        using namespace nonlocal::mesh::utils;

        conductivity.compute(
            parameters, 
            utils::inner_nodes(mesh->container(), boundaries_conditions), 
            is_symmetric, is_neumann,
            is_sol_depend ? std::optional{nodes_to_qnodes(*mesh, temperature_prev)} : std::nullopt
        );

        convection_condition_2d(conductivity.matrix_inner(), *mesh, boundaries_conditions);
        boundary_condition_first_kind_2d(f, *mesh, boundaries_conditions, conductivity.matrix_bound());

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

        if (!is_sol_depend)
            break;

        ++iteration;

        difference = (temperature_curr - temperature_prev).norm() / (temperature_curr.norm() ?: T{1});
        std::cout << "norm(prev - curr) = " << difference << std::endl;
    }
    std::cout << "Iterations: " << iteration << std::endl;
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Time: " << elapsed_seconds.count() << 's' << std::endl;

    return heat_equation_solution_2d<T, I>{mesh, parameters, temperature_curr};
}

}

#endif
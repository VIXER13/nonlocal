#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "boundary_condition_1d.hpp"
#include "heat_equation_parameters_1d.hpp"
#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"

//#include "thermal_conductivity_matrix_1d.hpp"
//#include "convection_condition_1d.hpp"
//#include <iostream>

namespace nonlocal::thermal {

template<class T>
constexpr bool is_neumann_problem(const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                  const std::array<T, 2>& alpha) noexcept {
    const bool left_neumann  = boundary_condition.front().type == boundary_condition_t::FLUX ||
                               boundary_condition.front().type == boundary_condition_t::CONVECTION &&
                               std::abs(alpha.front()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
    const bool right_neumann = boundary_condition.back().type == boundary_condition_t::FLUX ||
                               boundary_condition.back().type == boundary_condition_t::CONVECTION &&
                               std::abs(alpha.back()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
    return left_neumann && right_neumann;
}

template<class T>
constexpr bool is_robin_problem(const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition) noexcept {
    return boundary_condition.front().type == boundary_condition_t::CONVECTION &&
           boundary_condition.front().type == boundary_condition_t::CONVECTION;
}

template<class T>
constexpr bool is_solvable_robin_problem(const std::array<T, 2>& section, const std::array<T, 2>& alpha) noexcept {
    const T length = section.back() - section.front();
    return std::abs(alpha.back() / (length * alpha.back() - T{1}) - alpha.front()) > ROBIN_PROBLEM_ALPHA_THRESHOLD<T>;
}

template<class T, class I, class Right_Part>
std::vector<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                  const std::vector<equation_parameters<1, T, stationary_equation_parameters_1d>>& parameters,
                                                  const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                                  const Right_Part& right_part, const T energy = T{0}) {
    const bool is_neumann = false;
    //const bool is_neumann = is_neumann_problem(boundary_condition, equation_param.alpha);
    //if (is_neumann && std::abs(boundary_condition.front().val + boundary_condition.back().val) > NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>)
    //    throw std::domain_error{"Unsolvable Neumann problem: left_flow + right_flow != 0."};

    //const std::array<T, 2>& alpha = equation_param.alpha;
    //if (is_robin_problem(boundary_condition) && !is_solvable_robin_problem(mesh->section(), equation_param.alpha))
    //    throw std::domain_error{"Unsolvable Robin problem: alpha[0] == alpha[1] / ((b - a) * alpha[1] - 1)."};

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    //integrate_right_part(f, *mesh, right_part);
    boundary_condition_second_kind_1d(f, boundary_condition, std::array{size_t{0}, size_t(f.size() - 1 - is_neumann)});
    //convection_condition_right_part_1d(f, boundary_condition, equation_param.alpha);
    
    double time = omp_get_wtime();
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(
        parameters,
        { boundary_condition.front().type == boundary_condition_t::TEMPERATURE, 
          boundary_condition.back ().type == boundary_condition_t::TEMPERATURE },
        is_neumann
    );
    //convection_condition_matrix_part_1d(conductivity.matrix_inner(), boundary_type(boundary_condition), equation_param.alpha);
    std::cout << "Matrix calc: " << omp_get_wtime() - time << std::endl;

    boundary_condition_first_kind_1d(f, boundary_condition, conductivity.matrix_bound());
    //if (is_neumann)
    //    f[f.size()-1] = equation_param.integral;

    time = omp_get_wtime();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.matrix_inner()};
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    std::cout << "SLAE: " << omp_get_wtime() - time << std::endl;
    std::cout << "iterations: " << solver.iterations() << std::endl;
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh->nodes_count())};
}

/*
template<class T, class I, class Right_Part, class Influence_Function>
std::vector<T> stationary_heat_equation_solver_1d(const equation_parameters_1d<T>& equation_param,
                                                  const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                  const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
                                                  const Right_Part& right_part,
                                                  const T p1,
                                                  const Influence_Function& influence_function,
                                                  const Eigen::Matrix<T, Eigen::Dynamic, 1>& x0 = {}) {
    const bool is_neumann = is_neumann_problem(boundary_condition, equation_param.alpha);
    if (is_neumann && std::abs(boundary_condition.front().val + boundary_condition.back().val) > NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>)
        throw std::domain_error{"Unsolvable Neumann problem: left_flow + right_flow != 0."};

    const std::array<T, 2>& alpha = equation_param.alpha;
    if (is_robin_problem(boundary_condition) && !is_solvable_robin_problem(mesh->section(), equation_param.alpha))
        throw std::domain_error{"Unsolvable Robin problem: alpha[0] == alpha[1] / ((b - a) * alpha[1] - 1)."};

    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(equation_param.lambda, p1, influence_function, boundary_type(boundary_condition), is_neumann);
    convection_condition_matrix_part_1d(conductivity.matrix_inner(), boundary_type(boundary_condition), equation_param.alpha);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(conductivity.matrix_inner().cols());
    integrate_right_part(f, *mesh, right_part);
    boundary_condition_second_kind_1d(f, boundary_condition, std::array{size_t{0}, size_t(f.size() - 1 - is_neumann)});
    convection_condition_right_part_1d(f, boundary_condition, equation_param.alpha);
    boundary_condition_first_kind_1d(f, boundary_condition, conductivity.matrix_bound());
    if (is_neumann)
        f[f.size()-1] = equation_param.integral;

    const double time = omp_get_wtime();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.matrix_inner()};
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    std::cout << "SLAE: " << omp_get_wtime() - time << std::endl;
    std::cout << "iterations: " << solver.iterations() << std::endl;
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh->nodes_count())};
}
*/

}

#endif
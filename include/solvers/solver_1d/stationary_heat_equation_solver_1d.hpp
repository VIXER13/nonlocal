#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "parameters_1d.hpp"
#include "right_part_1d.hpp"
#include "thermal_conductivity_matrix_1d.hpp"
#include "boundary_condition_third_kind_1d.hpp"

namespace nonlocal::thermal {

template<class T>
constexpr bool is_neumann_problem(const std::array<stationary_boundary_1d_t<T>, 2>& boundary_condition,
                                  const std::array<T, 2>& alpha) noexcept {
    const bool left_neumann  = boundary_type(boundary_condition.front()) == boundary_condition_t::SECOND_KIND ||
                               boundary_type(boundary_condition.front()) == boundary_condition_t::THIRD_KIND &&
                               std::abs(alpha.front()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
    const bool right_neumann = boundary_type(boundary_condition.back()) == boundary_condition_t::SECOND_KIND ||
                               boundary_type(boundary_condition.back()) == boundary_condition_t::THIRD_KIND &&
                               std::abs(alpha.back()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
    return left_neumann && right_neumann;
}

template<class T>
constexpr bool is_robin_problem(const std::array<stationary_boundary_1d_t<T>, 2>& boundary_condition) noexcept {
    return boundary_type(boundary_condition.front()) == boundary_condition_t::THIRD_KIND &&
           boundary_type(boundary_condition.front()) == boundary_condition_t::THIRD_KIND;
}

template<class T>
constexpr bool is_solvable_robin_problem(const std::array<T, 2>& section, const std::array<T, 2>& alpha) noexcept {
    const T length = section.back() - section.front();
    return std::abs(alpha.back() / (length * alpha.back() - T{1}) - alpha.front()) > ROBIN_PROBLEM_ALPHA_THRESHOLD<T>;
}

template<class T, class I, class Right_Part, class Influence_Function>
std::vector<T> stationary_heat_equation_solver_1d(const nonlocal_parameters_1d<T>& nonloc_param,
                                                  const heat_equation_parameters_1d<T>& equation_param,
                                                  const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                  const std::array<stationary_boundary_1d_t<T>, 2>& boundary_condition,
                                                  const Right_Part& right_part,
                                                  const Influence_Function& influence_function) {
    const bool is_neumann = is_neumann_problem(boundary_condition, equation_param.alpha);
    if (is_neumann && std::abs(boundary_value(boundary_condition.front()) + boundary_value(boundary_condition.back())) > NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>)
        throw std::domain_error{"This is unsolvable Neumann problem: left_flow + right_flow != 0"};

    const std::array<T, 2>& alpha = equation_param.alpha;
    if (is_robin_problem(boundary_condition) && !is_solvable_robin_problem(mesh->section(), equation_param.alpha))
        throw std::domain_error{"This is unsolvable Robin problem: alpha[0] == alpha[1] / ((b - a) * alpha[1] - 1)."};

    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(equation_param.lambda, nonloc_param.p1, influence_function, boundary_type(boundary_condition), is_neumann);
    boundary_condition_third_kind_1d(conductivity.matrix_inner(), boundary_condition, equation_param.alpha);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(conductivity.matrix_inner().cols());
    integrate_right_part(f, *mesh, right_part);
    boundary_condition_second_kind_1d(f, boundary_condition, std::array<size_t, 2>{0, size_t(f.size() - 1 - is_neumann)});
    boundary_condition_first_kind_1d(f, boundary_condition, conductivity.matrix_bound());
    if (is_neumann)
        f[f.size()-1] = equation_param.integral;

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{conductivity.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh->nodes_count())};
}

}

#endif
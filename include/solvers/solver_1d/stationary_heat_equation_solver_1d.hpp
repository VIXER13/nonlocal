#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_third_kind.hpp"
#include "parameters_1d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Right_Part, class Influence_Function>
std::vector<T> stationary_heat_equation_solver_1d(const nonlocal_parameters_1d<T>& nonloc_param,
                                                  const heat_equation_parameters_1d<T>& equation_param,
                                                  const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                  const std::array<std::pair<boundary_condition_t, T>, 2>& boundary_condition,
                                                  const Right_Part& right_part,
                                                  const Influence_Function& influence_function) {
    const bool is_neumann = boundary_condition.front().first == boundary_condition_t::SECOND_KIND &&
                            boundary_condition.back().first  == boundary_condition_t::SECOND_KIND;
    if (is_neumann && boundary_condition.front().second + boundary_condition.back().second > MAX_NEUMANN_PROBLEM_ERROR<T>)
        throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};

    thermal_conductivity_matrix_1d<T, I> matrix{mesh};
    matrix.template calc_matrix(equation_param.lambda, {boundary_condition.front().first, boundary_condition.back().first},
                                nonloc_param.p1, influence_function);
    boundary_condition_third_kind_1d(matrix.matrix_inner(), boundary_condition, equation_param.alpha);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(matrix.matrix_inner().cols());
    integrate_right_part(f, *mesh, right_part);
    boundary_condition_second_kind_1d(f, boundary_condition, std::array<size_t, 2>{0, size_t(f.size() - 1 - is_neumann)});
    boundary_condition_first_kind_1d(f, boundary_condition, matrix.matrix_bound());

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{matrix.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh->nodes_count())};
}

}

#endif
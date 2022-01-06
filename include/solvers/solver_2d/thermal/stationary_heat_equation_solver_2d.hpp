#ifndef STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include "parameters_2d.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "right_part_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index, material_t Material, class Right_Part, class Influence_Function>
std::vector<T> stationary_heat_equation_solver_2d(const equation_parameters<T, Material>& equation_param,
                                                  const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                  const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
                                                  const Right_Part& right_part,
                                                  const T p1,
                                                  const Influence_Function& influence_function) {
    const bool is_neumann = std::all_of(boundary_condition.cbegin(), boundary_condition.cend(),
        [](const auto& bound) noexcept { return bound.second.cond.type == boundary_condition_t::FLUX; });

    const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), boundary_condition);
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh_proxy};
    conductivity.template calc_matrix(equation_param, is_inner, p1, influence_function, is_neumann);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(conductivity.matrix_inner().cols());
    const std::array<size_t, 2> first_last_nodes = {mesh_proxy->first_node(), mesh_proxy->last_node()};
    integrate_right_part<1>(f, *mesh_proxy, right_part);
    boundary_condition_second_kind_2d(f, *mesh_proxy, first_last_nodes, boundary_condition);
    boundary_condition_first_kind_2d(f, mesh_proxy->mesh(), first_last_nodes, boundary_condition, conductivity.matrix_bound());
    if (is_neumann)
        f[f.size()-1] = equation_param.integral;

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{conductivity.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh_proxy->mesh().nodes_count())};
}

}

#endif
#ifndef STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "convection_condition_2d.hpp"
#include "right_part_2d.hpp"
#include "heat_equation_solution_2d.hpp"

namespace nonlocal::thermal {

template<class T, class I>
bool is_solvable_neumann_problem(const mesh::mesh_proxy<T, I>& mesh_proxy, const Eigen::Matrix<T, Eigen::Dynamic, 1>& f) {
    const T sum = std::accumulate(std::next(f.begin(), mesh_proxy.first_node()),
                                  std::next(f.begin(), mesh_proxy.last_node ()), T{0});
    return std::abs(MPI_utils::reduce(sum)) < NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>;
}

template<class T, class I, class Matrix_Index, material_t Material, class Right_Part, class Influence_Function>
solution<T, I> stationary_heat_equation_solver_2d(const equation_parameters_2d<T, Material>& equation_param,
                                                  const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                  const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
                                                  const Right_Part& right_part,
                                                  const T p1,
                                                  const Influence_Function& influence_function) {
    const auto bounds_types = boundary_type(boundary_condition);
    const bool is_neumann = std::all_of(bounds_types.cbegin(), bounds_types.cend(),
        [](const auto& bound) constexpr noexcept { return bound.second.front() == boundary_condition_t::FLUX; });
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh_proxy->mesh().nodes_count() + is_neumann);
    boundary_condition_second_kind_2d(f, *mesh_proxy, boundary_condition);
    if (is_neumann) {
        if (!is_solvable_neumann_problem(*mesh_proxy, f))
            throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
        f[f.size()-1] = equation_param.integral;
    }

    const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), bounds_types);
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh_proxy};
    conductivity.template calc_matrix<Material>(equation_param.lambda, is_inner, p1, influence_function, is_neumann);
    convection_condition_2d(conductivity.matrix_inner(), *mesh_proxy, bounds_types, equation_param.alpha);
    integrate_right_part<1>(f, *mesh_proxy, right_part);
    boundary_condition_first_kind_2d(f, *mesh_proxy, boundary_condition, conductivity.matrix_bound());

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{conductivity.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return solution<T, I>{mesh_proxy, temperature};
}

}

#endif
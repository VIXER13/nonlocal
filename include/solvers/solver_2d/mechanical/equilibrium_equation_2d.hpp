#ifndef NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP
#define NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP

#include "temperature_condition_2d.hpp"
#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include "stiffness_matrix_2d.hpp"
#include "right_part_2d.hpp"
#include "mechanical_solution.hpp"

namespace nonlocal::mechanical {

template<class T, class I, class Matrix_Index, class Right_Part, class Influence_Function>
mechanical::solution<T, I> equilibrium_equation(const mechanical::equation_parameters<T>& parameters,
                                                const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 2>>& boundary_condition,
                                                const Right_Part& right_part,
                                                const T p1,
                                                const Influence_Function& influence_function) {
    const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), boundary_type(boundary_condition));
    stiffness_matrix<T, I, Matrix_Index> stiffness{mesh_proxy};
    stiffness.template calc_matrix(hooke_matrix(parameters), is_inner, p1, influence_function);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(stiffness.matrix_inner().cols());
    const std::array<size_t, 2> first_last_nodes = {2 * mesh_proxy->first_node(), 2 * mesh_proxy->last_node()};
    integrate_right_part<2>(f, *mesh_proxy, right_part);
    boundary_condition_second_kind_2d(f, *mesh_proxy, first_last_nodes, boundary_condition);
    if (parameters.thermoelasticity)
        temperature_condition(f, *mesh_proxy, parameters, p1, influence_function);
    boundary_condition_first_kind_2d(f, mesh_proxy->mesh(), first_last_nodes, boundary_condition, stiffness.matrix_bound());

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{stiffness.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> displacement = solver.solve(f);
    return solution<T, I>{mesh_proxy, parameters, influence_function, displacement};
}

}

#endif
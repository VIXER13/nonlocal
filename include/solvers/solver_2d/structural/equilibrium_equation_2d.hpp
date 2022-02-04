#ifndef NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP
#define NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP

#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include "stiffness_matrix_2d.hpp"
#include "right_part_2d.hpp"
#include "mechanical_solution.hpp"

namespace nonlocal::mechanical {

/*
template<class T, class I, class Influence_Function>
void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                           const mesh::mesh_proxy<T, I>& mesh_proxy,
                           const mechanical::equation_parameters<T>& parameters,
                           const T p1,
                           const Influence_Function& influence_fun) {
    const T nu = parameters.poisson(),
            E  = parameters.young();
    const T factor = 0.5 * parameters.alpha * E / (1 - nu);
    const std::array<std::vector<T>, 2> gradient = mesh_proxy.template calc_gradient(parameters.delta_temperature);
    using namespace metamath::function;
    const std::array<std::vector<T>, 2> eps_T = { factor * mesh_proxy.approx_in_quad(gradient[0]),
                                                  factor * mesh_proxy.approx_in_quad(gradient[1]) };

#pragma omp parallel for default(none) shared(f, eps_T, p1, mesh_proxy)
    for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node)
        for(const I e : mesh_proxy.nodes_elements_map(node)) {
            const size_t i = mesh_proxy().global_to_local_numbering(e).find(node)->second;
            const std::array<T, 2> integral = p1 * integrate_temperature_loc(eps_T, e, i);
            using enum axis;
            f[2 * node + X] += integral[X];
            f[2 * node + Y] += integral[Y];
        }

    if (parameters.p1 < MAX_NONLOCAL_WEIGHT<T>) {
        const size_t p2 = 1 - parameters.p1;
#pragma omp parallel for default(none) shared(f, eps_T, p2, mesh_proxy) firstprivate(influence_fun)
        for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node)
            for(const I eL : mesh_proxy.nodes_elements_map(node)) {
                const size_t iL = mesh_proxy.global_to_local_numbering(eL).find(node)->second;
                for(const I eNL : mesh_proxy().neighbors(eL)) {
                    const std::array<T, 2> integral = integrate_temperature_nonloc(eps_T, eL, eNL, iL, influence_fun);
                    using enum axis;
                    f[2 * node + X] += p2 * integral[X];
                    f[2 * node + Y] += p2 * integral[Y];
                }
            }
    }
}
*/

template<class T, class I, class Matrix_Index, class Right_Part, class Influence_Function>
mechanical::solution<T, I> equilibrium_equation(const mechanical::equation_parameters<T>& parameters,
                                                const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 2>>& boundary_condition,
                                                const Right_Part& right_part,
                                                const T p1,
                                                const Influence_Function& influence_function) {
    const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), boundary_condition);
    stiffness_matrix<T, I, Matrix_Index> stiffness{mesh_proxy};
    stiffness.template calc_matrix(hooke_matrix(parameters), is_inner, p1, influence_function);

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(stiffness.matrix_inner().cols());
    const std::array<size_t, 2> first_last_nodes = {2 * mesh_proxy->first_node(), 2 * mesh_proxy->last_node()};
    integrate_right_part<2>(f, *mesh_proxy, right_part);
    boundary_condition_second_kind_2d(f, *mesh_proxy, first_last_nodes, boundary_condition);
    boundary_condition_first_kind_2d(f, mesh_proxy->mesh(), first_last_nodes, boundary_condition, stiffness.matrix_bound());

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{stiffness.matrix_inner()};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> displacement = solver.solve(f);
    return solution<T, I>{mesh_proxy, parameters, influence_function, displacement};
}

}

#endif
#ifndef THERMOMECHANICAL_EQUATION_SOLVER_HPP
#define THERMOMECHANICAL_EQUATION_SOLVER_HPP

#include "static_equation_solver.hpp"

namespace thermomechanical_equation_with_nonloc
{

template<class Type>
struct parameters
{
    Type nu    = 0, // Коэффициент Пуассона
         E     = 0, // Модуль Юнга
         alpha = 0; // Коэффициент теплового линейного расширения
};

template<class Type, class Index, class Vector>
static void temperature_condition(const mesh_2d<Type, Index>& mesh, const parameters<Type>& params,
                                  const Vector& temperature_field, Vector& f)
{
    const Type coeff = params.alpha * params.E / (1 - 2 * params.nu);
    matrix<Type> coords, jacobi_matrices;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        approx_quad_nodes_coords(mesh, e, el, coords);
        approx_jacobi_matrices(mesh, e, el, jacobi_matrices);
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            Type integral_x = 0, integral_y = 0;
            for(size_t q = 0; q < e->qnodes_count(); ++q)
            {
                Type Txi = 0, Teta = 0,
                     jacobian = e->weight(q) * e->qN(i, q) * (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2));
                for(size_t j = 0; j < e->nodes_count(); ++j)
                {
                    Txi  += temperature_field[mesh.node_number(el, j)] * e->qNxi(j, q);
                    Teta += temperature_field[mesh.node_number(el, j)] * e->qNeta(j, q);
                }
                integral_x += Txi  * jacobian;
                integral_y += Teta * jacobian;
            }
            f[2*mesh.node_number(el, i)]   += coeff * integral_x;
            f[2*mesh.node_number(el, i)+1] += coeff * integral_y;
        }
    }
}

template<class Type, class Index>
static Eigen::Matrix<Type, Eigen::Dynamic, 1> 
    stationary(const mesh_2d<Type, Index>& mesh, const parameters<Type>& params,
               const std::vector<static_equation_with_nonloc::boundary_condition<Type>>& bounds_cond,
               const static_equation_with_nonloc::distributed_load<Type>& right_part,
               const Eigen::Matrix<Type, Eigen::Dynamic, 1>& temperature_field,
               const Type p1, const std::function<Type(Type, Type, Type, Type)>& influence_fun)
{
    using namespace static_equation_with_nonloc;

    Eigen::Matrix<Type, Eigen::Dynamic, 1> f = Eigen::Matrix<Type, Eigen::Dynamic, 1>::Zero(2*mesh.nodes_count());
    Eigen::SparseMatrix<Type, Eigen::ColMajor, Index> K      (2*mesh.nodes_count(), 2*mesh.nodes_count()),
                                                      K_bound(2*mesh.nodes_count(), 2*mesh.nodes_count());
    
    double time = omp_get_wtime();
    create_matrix(mesh, {.nu = params.nu, .E = params.E}, bounds_cond, K, K_bound, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    integrate_right_part(mesh, right_part, f);
    temperature_condition(mesh, params, temperature_field, f);

    time = omp_get_wtime();
    boundary_condition_calc(mesh, kinematic_nodes_vectors(mesh, bounds_cond), bounds_cond, K_bound, f);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::PardisoLDLT<Eigen::SparseMatrix<Type, Eigen::ColMajor, Index>, Eigen::Lower> solver;
    solver.compute(K);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;

    return u;
}

}

#endif
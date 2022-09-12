#ifndef NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP
#define NONLOCAL_EQUILIBRIUM_EQUATION_2D_HPP

#include "temperature_condition_2d.hpp"
#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include "stiffness_matrix_2d.hpp"
#include "right_part_2d.hpp"
#include "mechanical_solution_2d.hpp"
#include "conjugate_gradient.hpp"

namespace nonlocal::mechanical {

template<class T, class I, class Matrix_Index, class Right_Part, class Influence_Function>
mechanical::mechanical_solution_2d<T, I> equilibrium_equation(const mechanical::equation_parameters<T>& parameters,
                                                const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
                                                const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 2>>& boundary_condition,
                                                const Right_Part& right_part,
                                                const T p1,
                                                const Influence_Function& influence_function) {
    const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), boundary_type(boundary_condition));

    double time = omp_get_wtime();
    stiffness_matrix<T, I, Matrix_Index> stiffness{mesh_proxy};
    stiffness.template calc_matrix(hooke_matrix(parameters), is_inner, p1, influence_function);
    std::cout << "conduction matrix create: " << omp_get_wtime() - time << std::endl;
    std::cout << "inner matrix non-zero elements count: " << stiffness.matrix_inner().nonZeros() << std::endl;
    std::cout << "bound matrix non-zero elements count: " << stiffness.matrix_bound().nonZeros() << std::endl;

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(stiffness.matrix_inner().cols());
    integrate_right_part<2>(f, *mesh_proxy, right_part);
    boundary_condition_second_kind_2d(f, *mesh_proxy, boundary_condition);
    if (parameters.is_thermoelasticity)
        temperature_condition(f, *mesh_proxy, parameters, p1, influence_function);
    boundary_condition_first_kind_2d(f, *mesh_proxy, boundary_condition, stiffness.matrix_bound());

    time = omp_get_wtime();
    //const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{stiffness.matrix_inner()};
    //const Eigen::Matrix<T, Eigen::Dynamic, 1> displacement = solver.solve(f);
    const slae::conjugate_gradient<T, Matrix_Index> solver{stiffness.matrix_inner()};
    const auto displacement = solver.solve(f);
    std::cout << "Slae solve time: " << omp_get_wtime() - time << std::endl;
    std::cout << "iterations: " << solver.iterations() << std::endl;
    return mechanical_solution_2d<T, I>{mesh_proxy, p1, influence_function, parameters, displacement};
}

}

#endif
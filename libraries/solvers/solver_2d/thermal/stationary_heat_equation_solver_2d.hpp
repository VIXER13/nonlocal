#ifndef NONLOCAL_STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP
#define NONLOCAL_STATIONARY_HEAT_EQUATION_SOLVER_2D_HPP

#include "solvers_utils.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "thermal_boundary_conditions_2d.hpp"
#include "boundary_condition_first_kind_2d.hpp"
//#include "convection_condition_2d.hpp"
//#include "right_part_2d.hpp"
#include "heat_equation_solution_2d.hpp"

#include "conjugate_gradient.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class Matrix_Index, class T, class I, class Right_Part, class Influence_Function>
heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh,
                                                                   const nonlocal::thermal::parameter_2d<T>& parameters,
                                                                   const boundaries_conditions_2d<T>& boundaries_conditions,
                                                                   const Right_Part& right_part,
                                                                   const T p1,
                                                                   const Influence_Function& influence_function, 
                                                                   const T energy = T{0}) {
    static constexpr size_t DoF = 1;
    static constexpr auto is_second_kind = [](const boundary_condition_2d<T>& condition) {
        return bool(dynamic_cast<const flux_2d<T>*>(condition.get()));
    };
    const auto conditions = boundaries_conditions | std::views::values;
    const bool is_neumann = std::all_of(conditions.begin(), conditions.end(), is_second_kind);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->container().nodes_count() + is_neumann);
    // boundary_condition_second_kind_2d(f, *mesh_proxy, boundary_condition);
    // if (is_neumann) {
    //     if (!is_solvable_neumann_problem(*mesh_proxy, f))
    //         throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
    //     f[f.size()-1] = equation_param.integral;
    // }

    auto start_time = std::chrono::high_resolution_clock::now();
    thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh};
    conductivity.template compute(
        parameters.conductivity, parameters.material, 
        utils::inner_nodes<DoF>(mesh->container(), boundaries_conditions), 
        p1, influence_function, 
        is_neumann
    );
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Conductivity matrix calculated time: " << elapsed_seconds.count() << 's' << std::endl;
    // convection_condition_matrix_part_2d(conductivity.matrix_inner(), *mesh_proxy, bounds_types, equation_param.heat_transfer);
    // convection_condition_right_part_2d(f, *mesh_proxy, boundary_condition, equation_param.heat_transfer);
    // integrate_right_part<1>(f, *mesh_proxy, right_part);
    boundary_condition_first_kind_2d<DoF>(f, *mesh, boundaries_conditions, conductivity.matrix_bound());

    // std::cout << Eigen::MatrixXd{conductivity.matrix_inner()} << std::endl << std::endl;
    // std::cout << Eigen::MatrixXd{conductivity.matrix_bound()} << std::endl << std::endl;
    // std::cout << f.transpose() << std::endl << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    // const slae::conjugate_gradient<T, Matrix_Index> solver{conductivity.matrix_inner()};
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{conductivity.matrix_inner()};
    const auto temperature = solver.solve(f);
    // const auto temperature = f;
    // std::cout << temperature.transpose() << std::endl << std::endl;

    elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "SLAE solution time: " << elapsed_seconds.count() << 's' << std::endl;
    std::cout << "Iterations: " << solver.iterations() << std::endl;
    return heat_equation_solution_2d<T, I>{mesh, p1, influence_function, parameters, temperature};
}

// template<class T, class I>
// bool is_solvable_neumann_problem(const mesh::mesh_proxy<T, I>& mesh_proxy, const Eigen::Matrix<T, Eigen::Dynamic, 1>& f) {
//     const T sum = std::accumulate(std::next(f.begin(), mesh_proxy.first_node()),
//                                   std::next(f.begin(), mesh_proxy.last_node ()), T{0});
//     return std::abs(parallel_utils::reduce(sum)) < NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>;
// }

// template<class T, class I, class Matrix_Index, material_t Material, class Right_Part, class Influence_Function>
// heat_equation_solution_2d<T, I> stationary_heat_equation_solver_2d(const equation_parameters_2d<T, Material>& equation_param,
//                                                                    const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy,
//                                                                    const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
//                                                                    const Right_Part& right_part,
//                                                                    const T p1,
//                                                                    const Influence_Function& influence_function) {
//     const auto bounds_types = boundary_type(boundary_condition);
//     const bool is_neumann = std::all_of(bounds_types.cbegin(), bounds_types.cend(),
//         [](const auto& bound) constexpr noexcept { return bound.second.front() == boundary_condition_t::FLUX; });
//     Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh_proxy->mesh().nodes_count() + is_neumann);
//     boundary_condition_second_kind_2d(f, *mesh_proxy, boundary_condition);
//     if (is_neumann) {
//         if (!is_solvable_neumann_problem(*mesh_proxy, f))
//             throw std::domain_error{"Unsolvable Neumann problem: contour integral != 0."};
//         f[f.size()-1] = equation_param.integral;
//     }

//     const std::vector<bool> is_inner = inner_nodes(mesh_proxy->mesh(), bounds_types);
//     double time = omp_get_wtime();
//     thermal_conductivity_matrix_2d<T, I, Matrix_Index> conductivity{mesh_proxy};
//     conductivity.template calc_matrix<Material>(equation_param.thermal_conductivity, is_inner, p1, influence_function, is_neumann);
//     std::cout << "conduction matrix create: " << omp_get_wtime() - time << std::endl;
//     std::cout << "inner matrix non-zero elements count: " << conductivity.matrix_inner().nonZeros() << std::endl;
//     std::cout << "bound matrix non-zero elements count: " << conductivity.matrix_bound().nonZeros() << std::endl;
//     convection_condition_matrix_part_2d(conductivity.matrix_inner(), *mesh_proxy, bounds_types, equation_param.heat_transfer);
//     convection_condition_right_part_2d(f, *mesh_proxy, boundary_condition, equation_param.heat_transfer);
//     integrate_right_part<1>(f, *mesh_proxy, right_part);
//     boundary_condition_first_kind_2d(f, *mesh_proxy, boundary_condition, conductivity.matrix_bound());

//     const slae::conjugate_gradient<T, Matrix_Index> solver{conductivity.matrix_inner()};
//     const auto temperature = solver.solve(f);
//     std::cout << "Slae solve time: " << omp_get_wtime() - time << std::endl;
//     std::cout << "iterations: " << solver.iterations() << std::endl;
//     return heat_equation_solution_2d<T, I>{mesh_proxy, p1, influence_function, equation_param, temperature};
// }

}

#endif
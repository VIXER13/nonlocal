#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "mesh_1d_utils.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class T, class I, class Right_Part>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const parameters_1d<T>& parameters,
                                                                const thermal_boundaries_conditions_1d<T>& boundaries_conditions,
                                                                const Right_Part& right_part,
                                                                const T energy = T{0}) {
    static constexpr auto is_flux = [](const auto& condition) {
        return bool(dynamic_cast<const flux_1d<T>*>(condition.get()));
    };
    const bool is_neumann = std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_flux);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(f, boundaries_conditions, is_neumann);
    if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(f, *mesh, right_part);
    if (is_neumann && std::abs(std::reduce(f.begin(), f.end())) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};

    auto start_time = std::chrono::high_resolution_clock::now();
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    conductivity.template calc_matrix(
        parameters,
        { bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.front().get())),
          bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.back ().get())) },
        is_neumann
    );
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Conductivity matrix calculated time: " << elapsed_seconds.count() << 's' << std::endl;

    start_time = std::chrono::high_resolution_clock::now();
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature;
    if (is_neumann) {
        f[f.size() - 1] = energy;
        const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
        std::cout << "Iterations: " << solver.iterations() << std::endl;
    } else {
        convection_condition_1d(conductivity.matrix_inner(), boundaries_conditions);
        boundary_condition_first_kind_1d(f, conductivity.matrix_bound(), boundaries_conditions);
        const Eigen::SimplicialCholesky<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
    }
    elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "SLAE solution time: " << elapsed_seconds.count() << 's' << std::endl;
    return heat_equation_solution_1d<T>{mesh, parameters, temperature};
}



template<class T, class I, class Right_Part, class Initial_distribution>
heat_equation_solution_1d<T> stationary_nonlinear_heat_equation_solver_1d  (const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                            const parameters_1d<T>& parameters,
                                                                            const thermal_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                                            const Right_Part& right_part,
                                                                            const Initial_distribution& initial_dist,
                                                                            const T energy = T{0}) {
    static constexpr auto is_flux = [](const auto& condition) {
        return bool(dynamic_cast<const flux_1d<T>*>(condition.get()));
    };
    const bool is_neumann = std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_flux);
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    boundary_condition_second_kind_1d(f, boundaries_conditions, is_neumann);
    if constexpr (!std::is_same_v<Right_Part, std::remove_cvref_t<decltype(EMPTY_FUNCTION)>>)
        integrate_right_part(f, *mesh, right_part);

    if (is_neumann && std::abs(std::reduce(f.begin(), f.end())) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};

    // Запуск итерационного процесса
    if (is_neumann) 
        f[f.size() - 1] = energy;

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count());
    for(size_t i = 0; i < temperature_prev.size(); ++i) 
        temperature_prev[i] = initial_dist;

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = temperature_prev;
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_f = f;
    size_t iters = 0;
    const size_t maxiters = temperature_curr.size() * 1.3;
    T norm = 0;

    using namespace nonlocal::mesh::utils;

    auto start_time = std::chrono::high_resolution_clock::now();
    do{
        iters++;
        std::swap(temperature_prev, temperature_curr);

        //Проецируем решение с nodes на qnodes
        // std::vector<T> sol_in_qnodes(mesh->nodes_count(), T{0});
        // std::copy(temperature_prev.begin(), temperature_prev.end(), sol_in_qnodes.begin());
        const std::vector<T> sol_in_qnodes = from_nodes_to_qnodes(*mesh, temperature_prev);

        // Откатываем правую часть к изначальному значению
        std::copy(initial_f.begin(), initial_f.end(), f.begin());

        // На каждой итерации корректируем (пересобираем) матрицу, хотим сойтись к точному решению стационарной задачи
        conductivity.template calc_matrix(
            parameters,
            { bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.front().get())),
              bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.back ().get())) },
            is_neumann,
            sol_in_qnodes
        );

        if (is_neumann) {
            const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.matrix_inner()};
            temperature_curr = solver.solve(f);
        } else {
            convection_condition_1d(conductivity.matrix_inner(), boundaries_conditions);
            boundary_condition_first_kind_1d(f, conductivity.matrix_bound(), boundaries_conditions);
            const Eigen::SimplicialCholesky<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>> solver{conductivity.matrix_inner()};
            temperature_curr = solver.solve(f);
        }

        norm = (temperature_curr - temperature_prev).norm();
        std::cout << "norm(prev - curr) = " << norm << std::endl;
    } while ((iters < maxiters) && (norm > 1e-10));
    std::cout << "Iterations: " << iters << std::endl;
    std::chrono::duration<double> elapsed_seconds = std::chrono::high_resolution_clock::now() - start_time;
    std::cout << "Time: " << elapsed_seconds.count() << 's' << std::endl;
    std::cout << "Mean iteration time: " << elapsed_seconds.count() / iters << 's' << std::endl;

    return heat_equation_solution_1d<T>{mesh, parameters, temperature_curr};
   
}

}

#endif
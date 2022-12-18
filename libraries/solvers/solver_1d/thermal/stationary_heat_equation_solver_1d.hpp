#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

#include <chrono>

namespace nonlocal::thermal {

template<class T, class I, class Right_Part>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const std::vector<equation_parameters<1, T, parameters_1d>>& parameters,
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

}

#endif
#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "convection_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

namespace nonlocal::thermal {

template<class T>
constexpr bool is_neumann_problem(const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition) noexcept {
    return dynamic_cast<stationary_flux_1d<T>*>(boundary_condition.front().get()) &&
           dynamic_cast<stationary_flux_1d<T>*>(boundary_condition.back ().get());
}

template<class T, class I, class Right_Part>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const std::vector<equation_parameters<1, T, parameters_1d>>& parameters,
                                                                const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                                                                const Right_Part& right_part, 
                                                                const T energy = T{0}) {
    const bool is_neumann = is_neumann_problem<T>(boundary_condition);
    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    double time = omp_get_wtime();
    conductivity.template calc_matrix(
        parameters,
        { bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.front().get())),
          bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.back ().get())) },
        is_neumann
    );
    std::cout << "Matrix time: " << omp_get_wtime() - time << std::endl;

    const std::array<size_t, 2> indices = {0, mesh->nodes_count() - 1};
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    integrate_right_part(f, *mesh, right_part);
    for(const size_t b : std::ranges::iota_view{0u, boundary_condition.size()}) {
        convection_condition_1d(conductivity.matrix_inner(), *boundary_condition[b], indices[b]);
        boundary_condition_second_kind_1d<T>(f, *boundary_condition[b], indices[b]);
        boundary_condition_first_kind_1d<T>(f, conductivity.matrix_bound()[b], *boundary_condition[b], indices[b]);
    }
    if (is_neumann)
        f[f.size() - 1] = energy;

    time = omp_get_wtime();
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature;
    if (is_neumann) {
        const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
    } else {
        const Eigen::SimplicialCholesky<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>> solver{conductivity.matrix_inner()};
        temperature = solver.solve(f);
    }
    std::cout << "SLAE time: " << omp_get_wtime() - time << std::endl;
    return heat_equation_solution_1d<T>{mesh, parameters, temperature};
}

}

#endif
#ifndef STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP
#define STATIONARY_HEAT_EQUATION_SOLVER_1D_HPP

#include "heat_equation_parameters_1d.hpp"
#include "thermal_conductivity_matrix_1d.hpp"
#include "right_part_1d.hpp"
#include "boundary_condition_first_kind_1d.hpp"
#include "boundary_condition_second_kind_1d.hpp"
#include "heat_equation_solution_1d.hpp"
#include "thermal_boundary_conditions_1d.hpp"
#include "convection_condition_1d.hpp"

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_1d final {
    T conductivity = T{1};
};

// template<class T>
// constexpr bool is_neumann_problem(const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition,
//                                   const std::array<T, 2>& alpha) noexcept {
//     const bool left_neumann  = boundary_condition.front().type == boundary_condition_t::FLUX ||
//                                boundary_condition.front().type == boundary_condition_t::CONVECTION &&
//                                std::abs(alpha.front()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
//     const bool right_neumann = boundary_condition.back().type == boundary_condition_t::FLUX ||
//                                boundary_condition.back().type == boundary_condition_t::CONVECTION &&
//                                std::abs(alpha.back()) < NEUMANN_PROBLEM_ALPHA_THRESHOLD<T>;
//     return left_neumann && right_neumann;
// }

// template<class T>
// constexpr bool is_robin_problem(const std::array<stationary_boundary_1d_t<boundary_condition_t, T>, 2>& boundary_condition) noexcept {
//     return boundary_condition.front().type == boundary_condition_t::CONVECTION &&
//            boundary_condition.front().type == boundary_condition_t::CONVECTION;
// }

// template<class T>
// constexpr bool is_solvable_robin_problem(const std::array<T, 2>& section, const std::array<T, 2>& alpha) noexcept {
//     const T length = section.back() - section.front();
//     return std::abs(alpha.back() / (length * alpha.back() - T{1}) - alpha.front()) > ROBIN_PROBLEM_ALPHA_THRESHOLD<T>;
// }

template<class T>
constexpr bool is_neumann_problem(const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition) noexcept {
    return dynamic_cast<stationary_flux_1d<T>*>(boundary_condition.front().get()) &&
           dynamic_cast<stationary_flux_1d<T>*>(boundary_condition.back ().get());
}

template<class T, class I, class Right_Part>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const std::vector<equation_parameters<1, T, stationary_equation_parameters_1d>>& parameters,
                                                                const std::array<std::unique_ptr<stationary_thermal_boundary_condition_1d>, 2>& boundary_condition,
                                                                const Right_Part& right_part, 
                                                                const T energy = T{0}) {
    const bool is_neumann = is_neumann_problem<T>(boundary_condition);
    //if (is_neumann && std::abs(boundary_condition.front()() + boundary_condition.back()()) > NEUMANN_PROBLEM_MAX_BOUNDARY_ERROR<T>)
    //   throw std::domain_error{"Unsolvable Neumann problem: left_flow + right_flow != 0."};

    //const std::array<T, 2>& alpha = equation_param.alpha;
    //if (is_robin_problem(boundary_condition) && !is_solvable_robin_problem(mesh->section(), equation_param.alpha))
    //   throw std::domain_error{"Unsolvable Robin problem: alpha[0] == alpha[1] / ((b - a) * alpha[1] - 1)."};

    thermal_conductivity_matrix_1d<T, I> conductivity{mesh};
    double time = omp_get_wtime();
    conductivity.template calc_matrix(
        parameters,
        { bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.front().get())),
          bool(dynamic_cast<stationary_temperature_1d<T>*>(boundary_condition.back ().get())) },
        is_neumann
    );
    std::cout << "Matrix time: " << omp_get_wtime() - time << std::endl;

    const std::array<size_t, 2> indices = {0, size_t(conductivity.matrix_inner().cols() - 1)};
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + is_neumann);
    integrate_right_part(f, *mesh, right_part);
    for(const size_t b : std::ranges::iota_view{size_t{0}, boundary_condition.size()}) {
        convection_condition_1d(conductivity.matrix_inner(), *boundary_condition[b], indices[b]);
        boundary_condition_second_kind_1d<T>(f, *boundary_condition[b], indices[b]);
        boundary_condition_first_kind_1d<T>(f, conductivity.matrix_bound()[b], *boundary_condition[b], indices[b]);
    }
    if (is_neumann)
        f[f.size()-1] = energy;

    time = omp_get_wtime();
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>> solver{conductivity.matrix_inner()};
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    std::cout << "SLAE time: " << omp_get_wtime() - time << std::endl;
    return heat_equation_solution_1d<T>{mesh, parameters, temperature};
}

}

#endif
#pragma once

#include "thermal_conductivity_assembler.hpp"
#include "convection_condition_1d.hpp"
#include "radiation_condition_1d.hpp"
#include "heat_equation_solution_1d.hpp"

#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/base/right_part_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>

namespace nonlocal::thermal {

template<class T>
struct stationary_equation_parameters_1d final {
    std::optional<std::function<T(const T)>> right_part;
    std::optional<std::function<T(const T)>> initial_distribution;
    T tolerance = std::is_same_v<T, float> ? 1e-6 : 1e-15;
    size_t max_iterations = 100;
    T energy = T{0};
};

struct problem_settings final {
    bool is_neumann = false;
    bool is_nonlocal = false;
    bool is_nonconstant_parameters = false;
    bool is_radiation_boundary = false;
    bool is_solution_dependent = false;

    constexpr bool is_nonlinear() const noexcept {
        return is_radiation_boundary || is_solution_dependent;
    }

    constexpr bool is_symmetric() const noexcept {
        return !(is_nonlocal && is_nonconstant_parameters);
    }
};

template<class T>
problem_settings init_problem_settings(const parameters_1d<T>& parameters,
                                       const thermal_boundaries_conditions_1d<T>& boundaries_conditions) {
    static constexpr auto is_flux = [](const auto& condition) noexcept {
        return  bool(dynamic_cast<const flux_1d<T>*>(condition.get())) &&
               !bool(dynamic_cast<const combined_flux_1d<T>*>(condition.get()));
    };
    static constexpr auto is_nonlocal = [](const theory_t theory) noexcept {
        return theory == theory_t::NONLOCAL;
    };
    static constexpr auto is_radiation_boundary = [](const auto& condition) noexcept {
        return bool(dynamic_cast<const radiation_1d<T>*>(condition.get()));
    };
    static constexpr auto is_nonconstant_parameters = [](const auto& parameter) noexcept {
        return std::holds_alternative<spatial_dependency<T, 1>>(parameter.physical.conductivity) ||
               std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.conductivity);
    };
    static constexpr auto is_solution_dependent = [](const auto& parameter) noexcept {
        return std::holds_alternative<solution_dependency<T, 1>>(parameter.physical.conductivity);
    };

    const std::vector<theory_t> theories = theories_types(parameters);
    return {
        .is_neumann = std::all_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_flux),
        .is_nonlocal = std::any_of(theories.begin(), theories.end(), is_nonlocal),
        .is_nonconstant_parameters = std::any_of(parameters.begin(), parameters.end(), is_nonconstant_parameters),
        .is_radiation_boundary = std::any_of(boundaries_conditions.begin(), boundaries_conditions.end(), is_radiation_boundary),
        .is_solution_dependent = std::any_of(parameters.begin(), parameters.end(), is_solution_dependent)
    };
}

template<class T, class I>
heat_equation_solution_1d<T> stationary_heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                                                const parameters_1d<T>& parameters,
                                                                const thermal_boundaries_conditions_1d<T>& boundaries_conditions,  
                                                                const stationary_equation_parameters_1d<T>& additional_parameters) {
    const auto settings = init_problem_settings(parameters, boundaries_conditions);
    logger::info() << (settings.is_symmetric() ? "Symmetric problem" : "Asymmetrical problem") << std::endl;
    if (settings.is_neumann)
        logger::info() << "Neuman problem" << std::endl;
    if (settings.is_nonlinear())
        logger::info() << "Nonlinear problem" << std::endl;
    if (settings.is_solution_dependent)
        throw std::domain_error{"Parametrically nonlinear problems are not supported at the moment."};

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh->nodes_count() + settings.is_neumann);
    boundary_condition_second_kind_1d(f, boundaries_conditions, settings.is_neumann);
    if (additional_parameters.right_part)
        integrate_right_part(f, *mesh, *additional_parameters.right_part);
    if (settings.is_neumann && std::abs(std::reduce(f.begin(), f.end(), T{0})) > NEUMANN_PROBLEM_ERROR_THRESHOLD<T>)
        throw std::domain_error{"It's unsolvable Neumann problem."};
    if (settings.is_neumann) 
        f[f.size() - 1] = additional_parameters.energy;
    const Eigen::Matrix<T, Eigen::Dynamic, 1> initial_f = f;
    Eigen::Matrix<T, Eigen::Dynamic, 1> residual = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(f.size());

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_prev = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(f.size());
    if (settings.is_nonlinear() && additional_parameters.initial_distribution)
        for(const size_t node : mesh->nodes())
            temperature_prev[node] = (*additional_parameters.initial_distribution)(mesh->node_coord(node));
    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature_curr = temperature_prev;
    
    const std::array<bool, 2> is_first_kind = { 
        bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.front().get())),
        bool(dynamic_cast<temperature_1d<T>*>(boundaries_conditions.back ().get()))
    };
    finite_element_matrix_1d<T, I> conductivity;
    init_matrix_portrait(conductivity.inner, *mesh, theories_types(parameters), is_first_kind, settings.is_neumann, settings.is_symmetric());
    thermal_conductivity_assembler_1d<T, I> assembler{conductivity, mesh};

    T difference = T{1};
    size_t iteration = 0;
    do {
        if (settings.is_nonlinear()) {
            std::swap(temperature_prev, temperature_curr);
            std::copy(initial_f.begin(), initial_f.end(), f.begin());
        }
        using nonlocal::mesh::utils::from_nodes_to_qnodes;
        conductivity.set_zero();
        assembler.calc_matrix(
            parameters, is_first_kind, settings.is_neumann, settings.is_symmetric(),
            settings.is_solution_dependent ? std::optional{from_nodes_to_qnodes(*mesh, temperature_prev)} : std::nullopt
        );

        if (settings.is_neumann) {
            if (settings.is_symmetric()) {
                const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{conductivity.inner};
                temperature_curr = solver.solveWithGuess(f, temperature_prev);
            } else {
                const Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor, I>> solver{conductivity.inner};
                temperature_curr = solver.solveWithGuess(f, temperature_prev);
            }
        } else {
            convection_condition_1d(conductivity.inner, boundaries_conditions);
            boundary_condition_first_kind_1d(f, conductivity.bound, boundaries_conditions);
            if (settings.is_radiation_boundary) {
                if (settings.is_symmetric())
                    residual = conductivity.inner.template selfadjointView<Eigen::Upper>() * temperature_prev - f;
                else
                    residual = conductivity.inner * temperature_prev - f;
                radiation_condition_1d(conductivity.inner, residual, boundaries_conditions, temperature_prev);
            }

            if (settings.is_symmetric()) {
                const Eigen::SimplicialCholesky<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper, Eigen::NaturalOrdering<I>
                > solver{conductivity.inner};
                if (settings.is_radiation_boundary)
                    temperature_curr = temperature_prev - solver.solve(residual);
                else
                    temperature_curr = solver.solve(f);
            } else {
                const Eigen::SparseLU<
                    Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::NaturalOrdering<I>
                > solver{conductivity.inner};
                if (settings.is_radiation_boundary)
                    temperature_curr = temperature_prev - solver.solve(residual);
                else
                    temperature_curr = solver.solve(f);
            }
        }
        ++iteration;
        difference = (temperature_curr - temperature_prev).norm() / (temperature_curr.norm() ?: T{1});
    } while(settings.is_nonlinear() &&
            iteration < additional_parameters.max_iterations && 
            difference > additional_parameters.tolerance);
    return heat_equation_solution_1d<T>{mesh, parameters, temperature_curr};
}

}
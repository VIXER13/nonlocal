#pragma once

#include "save_results.hpp"

#include <config/task_data.hpp>
#include <config/mechanical_auxiliary_data.hpp>
#include <config/time_data.hpp>
#include <config/read_mesh.hpp>
#include <config/frequency_data.hpp>
#include <config/read_mechanical_boundary_conditions.hpp>
#include <config/read_mechanical_parameters.hpp>
#include <solvers/solver_1d/mechanical/stationary_mechanical_equation_solver_1d.hpp>
#include <solvers/solver_1d/mechanical/harmonic_mechanical_equation_solver_1d.hpp>

#include <chrono>

namespace nonlocal {

template<std::floating_point T>
void save_solution(const solver_1d::mechanical::mechanical_equation_solution_1d<T>& solution, 
                   const config::save_data& save,
                   const std::optional<uint64_t> step = std::nullopt) {
    if (step.has_value())
        logger::info() << "save step " << *step << std::endl;
    const std::filesystem::path path = step ? save.make_path(save.get_name("csv", "solution") + "_" + std::to_string(*step), "csv") : 
                                              save.path("csv", "csv", "solution");
    mesh::utils::save_as_csv(path, solution.mesh(), {{"displacement", solution.displacement()}, {"stress", solution.stress()}}, save.precision());
}

template<std::floating_point T, std::integral I>
void solve_mechanical_1d_problem(const nlohmann::json& config, const config::save_data& save, const config::analysis_type_t analysis_type) {
    const auto mesh = config::read_mesh_1d<T>(config, {});
    const auto parameters = config::read_mechanical_parameters_1d<T>(config["materials"], "materials");
    const auto auxiliary = config::mechanical_auxiliary_data_1d<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
    const auto boundaries_conditions = config::read_mechanical_boundaries_conditions_1d<T>(config["boundaries"], "boundaries");
    const auto right_part_input = [right_part = auxiliary.right_part](T x) {
        return std::visit(metamath::visitor{
            [](const T value) { return value; },
            [&x](const spatial_dependency<T, 1>& value) { return value(x); },
            [](const auto&) { throw std::domain_error{"Unsuported right part format."}; return T{0}; }
        }, right_part);
    };
    const auto initial_distribution_input = [value = auxiliary.initial_distribution](const T x) constexpr noexcept { return value; };

    using namespace solver_1d::mechanical;
    switch (analysis_type) {
        case config::analysis_type_t::Stationary: {
            auto solution = stationary_mechanical_equation_solver_1d<T, I>(
                mesh, parameters, boundaries_conditions,
                stationary_equation_parameters_1d<T>{
                    .right_part = right_part_input,
                    .initial_distribution = initial_distribution_input
                }
            );
            solution.calc_stress();
            save_solution(solution, save);
            break;
        }
        case config::analysis_type_t::TimeDependent: 
            throw std::domain_error{"TimeDependent analysis type for mechanical one-dimensional problem is not supported."};
        case config::analysis_type_t::TimeHarmonic: { 
            config::check_required_fields(config, {"frequency"});
            const config::frequency_data<T> sweep{config["frequency"], "frequency"};
            for (size_t i = 0; i < sweep.frequencies.size(); ++i) {
                auto start = std::chrono::steady_clock::now();
                auto solution = harmonic_mechanical_equation_solver_1d<T, I>(
                    mesh, parameters, boundaries_conditions,
                    time_harmonic_equation_parameters_1d<T>{
                        .right_part = right_part_input,
                        .initial_distribution = initial_distribution_input,
                        .frequency = sweep.frequencies[i]
                    }
                );
                solution.calc_stress();
                save_solution(solution, save, i);
                auto finish = std::chrono::steady_clock::now();
                logger::info() << " freq #"     << i << ": " << sweep.frequencies[i] << 
                                  " | time = " << std::chrono::duration_cast<std::chrono::seconds>(finish - start).count() << " sec." << std::endl;
            }   
            break;
        }
        case config::analysis_type_t::Unknown: 
        default: 
            throw std::domain_error{"Unknown analysis type for mechanical one-dimensional problem."};
    }
}

}
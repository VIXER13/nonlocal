#pragma once

#include "thermal_problems_1d.hpp"
#include "mechanical_problems_1d.hpp"
#include "save_results.hpp"

#include <config/read_mechanical_boundary_conditions.hpp>
#include <config/read_mechanical_parameters.hpp>
#include <config/read_mesh.hpp>
#include <config/read_thermal_boundary_conditions.hpp>
#include <config/task_data.hpp>
#include <config/save_data.hpp>
#include <config/time_data.hpp>
#include <config/thermal_auxiliary_data.hpp>
#include <mesh/mesh_2d/find_neighbours.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>
#include <solvers/solver_2d/thermal/nonstationary_heat_equation_solver_2d.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>
#include <solvers/solver_2d/mechanical/motion_equation_solver.hpp>

#include <set>

namespace nonlocal {

class _determine_problem final {
    static constexpr bool Only_Local = true;
    static constexpr bool Symmetric = true;

    constexpr explicit _determine_problem() noexcept = default;

    static void init_save_data(const config::save_data& save, const nlohmann::json& config);
    static std::vector<std::string> get_required_fields(const config::task_data& task);
    static bool is_thermal(const config::problem_t problem);
    static bool is_mechanical(const config::problem_t problem);

    template<std::floating_point T>
    friend void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T>
    friend void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T>
    friend void thermal_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::save_data& save);

    template<std::floating_point T>
    friend void mechanical_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::save_data& save);

    template<std::floating_point T>
    friend std::optional<solver_2d::thermal::heat_equation_solution_2d<T>> thermal_stationary_2d(
        std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem);

    template<std::floating_point T>
    friend std::optional<solver_2d::mechanical::mechanical_solution_2d<T>> mechanical_2d(
        std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem, const std::vector<T>& delta_temperature);

public:
    template<std::floating_point T>
    friend void determine_problem(const nlohmann::json& config);
};

template<std::floating_point T>
void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    if (task.problem == config::problem_t::Unknown)
        throw std::domain_error{"Unknown task. In the one-dimensional case, the following problems are available: "
                                "\"thermal\" and \"mechanical\""};
    if (parallel::MPI_rank() != 0) {
        logger::warning() << "Calculations are available only on the master process. "
                             "The current process has completed its work." << std::endl;
        return;
    }
    config::check_required_fields(config, {"boundaries", "materials"});
    config::check_optional_fields(config, {"mesh", "auxiliary"});
    switch (task.problem) {
        case config::problem_t::Thermal: {
            solve_thermal_1d_problem<T, int64_t>(config, save, task.analysis_type);
            break;
        }
        case config::problem_t::Mechanical: {
            solve_mechanical_1d_problem<T, int64_t>(config, save, task.analysis_type);
            break;
        }
        default: break;
    }
}

template<std::floating_point T>
std::optional<solver_2d::thermal::heat_equation_solution_2d<T>> thermal_stationary_2d(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem) {
    using DP = _determine_problem;
    if (!DP::is_thermal(problem))
        return std::nullopt;
    mesh->neighbours(mesh::find_neighbours(*mesh, config::read_influences<T>(config["materials"], "materials", "thermal")));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::Memory, !DP::Only_Local, DP::Symmetric);
    const auto boundaries_field = problem == config::problem_t::Thermal ? "boundaries" : "thermal_boundaries";
    return solver_2d::thermal::stationary_heat_equation_solver_2d(mesh,
        config::read_thermal_parameters_2d<T>(config["materials"], "materials"), 
        config::read_thermal_boundaries_conditions_2d<T>(config[boundaries_field], boundaries_field), 
        config::read_stationary_equation_parameters_2d<T>(config.value("auxiliary", nlohmann::json::object()), "auxiliary")
    );
}

template<std::floating_point T>
void thermal_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::save_data& save) {
    using DP = _determine_problem;
    mesh->neighbours(mesh::find_neighbours(*mesh, config::read_influences<T>(config["materials"], "materials", "thermal")));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::Memory, !DP::Only_Local, DP::Symmetric);
    const auto time = config::time_data<T>{config["time"], "time"};
    
    solver_2d::thermal::nonstationary_heat_equation_solver_2d<T> solver{mesh, time.time_step};
    const auto parameters = config::read_thermal_parameters_2d<T>(config["materials"], "materials");
    const auto boundaries_conditions = config::read_thermal_boundaries_conditions_2d<T>(config["boundaries"], "boundaries");
    const auto auxiliary = config::read_stationary_equation_parameters_2d<T>(config.value("auxiliary", nlohmann::json::object()), "auxiliary");
    const auto conductivity_parameters = evaluate_conductivity(*mesh, parameters, std::vector<T>(mesh->quad_shift(mesh->container().elements_2d_count()), T{0}));
    solver.compute(parameters, boundaries_conditions, auxiliary.initial_distribution);
    {
        solver_2d::thermal::heat_equation_solution_2d<T> solution{mesh, conductivity_parameters, solver.temperature()};
        solution.calc_flux();
        // save_solution(solution, save, 0u);
    }
    for(const uint64_t step : std::ranges::iota_view{1u, time.steps_count + 1}) {
        solver.calc_step(boundaries_conditions, auxiliary.right_part);
        if (step % time.save_frequency == 0) {
            logger::info() << "saving step " << step << std::endl;
            // solver_2d::thermal::heat_equation_solution_2d<T> solution{mesh, conductivity_parameters, solver.temperature()};
            // solution.calc_flux();
            // save_csv(solution, {}, save, step);
            // save_vtk(solution, {}, save, step);
        }
    }
}

template<std::floating_point T>
void mechanical_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::save_data& save) {
    using DP = _determine_problem;
    mesh->neighbours(mesh::find_neighbours(*mesh, config::read_influences<T>(config["materials"], "materials", "mechanical")));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::Memory, !DP::Only_Local, DP::Symmetric);
    constexpr auto Boundaries_Field = "boundaries";
    const config::time_data<T> time{config["time"], "time"};
    solver_2d::mechanical::motion_equation_solver<T> solver{mesh};
    solver.compute(config::read_mechanical_parameters_2d<T>(config["materials"], "materials"),
                   config::read_mechanical_boundaries_conditions_2d<T>(config[Boundaries_Field], Boundaries_Field),
                   time.time_step, time.initial_time);
    for(const size_t step : std::ranges::iota_view{0zu, time.steps_count}) {
        solver.calc_step();
        if (step % time.save_frequency == 0) {
            logger::info() << "saving step " << step << std::endl;
            const std::optional<solver_2d::mechanical::mechanical_solution_2d<T>> solution = solver.solution();
            save_csv({}, solution, save, step);
            save_vtk({}, solution, save, step);
        }
    }
}

template<std::floating_point T>
std::optional<solver_2d::mechanical::mechanical_solution_2d<T>> mechanical_2d(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem, const std::vector<T>& delta_temperature) {
    using DP = _determine_problem;
    if (!DP::is_mechanical(problem))
        return std::nullopt;
    mesh->neighbours(mesh::find_neighbours(*mesh, config::read_influences<T>(config["materials"], "materials", "mechanical")));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::Memory, !DP::Only_Local, DP::Symmetric);
    const auto boundaries_field = problem == config::problem_t::Mechanical ? "boundaries" : "mechanical_boundaries";
    return solver_2d::mechanical::equilibrium_equation(mesh, 
        config::read_mechanical_parameters_2d<T>(config["materials"], "materials"),
        config::read_mechanical_boundaries_conditions_2d<T>(config[boundaries_field], boundaries_field),
        delta_temperature
    );
}

template<std::floating_point T>
void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    if (task.problem == config::problem_t::Unknown)
        throw std::domain_error{"Unknown task. In the two-dimensional case, the following problems are available: "
                                "\"thermal\", \"mechanical\" and \"thermomechanical\""};
    
    using DP = _determine_problem;
    config::check_required_fields(config, DP::get_required_fields(task));
    config::check_optional_fields(config, {"auxiliary"});
    auto mesh = config::read_mesh_2d<T, uint32_t>(config["mesh"], "mesh");
    switch (task.analysis_type) {
        case config::analysis_type_t::Stationary: {
            const std::optional<solver_2d::thermal::heat_equation_solution_2d<T>> thermal_solution = thermal_stationary_2d<T>(mesh, config, task.problem);
            const std::optional<solver_2d::mechanical::mechanical_solution_2d<T>> mechanical_solution =
            mechanical_2d<T>(mesh, config, task.problem, thermal_solution ? thermal_solution->temperature() : std::vector<T>{});
            save_csv(thermal_solution, mechanical_solution, save);
            save_vtk(thermal_solution, mechanical_solution, save);
            break;
        }
        case config::analysis_type_t::Time_Harmonic: {
            throw std::domain_error{"Time_Harmonic analysis type for two-dimensional problem is not supported."};
        }
        case config::analysis_type_t::Time_Dependent: { 
            thermal_nonstationary_2d<T>(mesh, config, save);
            break;
        }
        case config::analysis_type_t::Unknown: 
        default: 
            throw std::domain_error{"Unknown analysis type for two-dimensional problem."};
    }
}

template<std::floating_point T>
void determine_problem(const nlohmann::json& config) {
    config::check_required_fields(config, {"task"});
    config::save_data save; 
    if (!config.contains("save"))
        logger::warning() << "There is no \"save\" field in the config. Required data may not be saved." << std::endl;
    else {
        save = config::save_data{config["save"], "save"};
        _determine_problem::init_save_data(save, config);
    }
    
    if (const config::task_data task{config["task"], "task"}; task.dimension == 1)
        problems_1d<T>(config, save, task);
    else if (task.dimension == 2)
        problems_2d<T>(config, save, task);
    else throw std::domain_error{"Problem dimension " + std::to_string(task.dimension) + " is not supported"};
}

}
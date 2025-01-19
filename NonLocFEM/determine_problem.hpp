#pragma once

#include "config/read_thermal_boundary_conditions.hpp"

#include "thermal_problems_1d.hpp"
#include "thermal_problems_2d.hpp"
#include "mechanical_problems_2d.hpp"

#include <set>

namespace nonlocal {

class _determine_problem final {
    static constexpr bool ONLY_LOCAL = true;
    static constexpr bool SYMMTERIC = true;

    constexpr explicit _determine_problem() noexcept = default;

    static void init_save_data(const config::save_data& save, const nlohmann::json& config);
    static std::vector<std::string> get_required_fields(const config::task_data& task);
    static bool is_thermal(const config::problem_t problem);
    static bool is_mechanical(const config::problem_t problem);

    template<std::floating_point T, std::signed_integral I>
    friend void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T, std::signed_integral I>
    friend void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task);

    template<std::floating_point T, std::signed_integral I>
    friend void thermal_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config,
                                         const config::save_data& save, const config::problem_t problem);

    template<std::floating_point T, std::signed_integral I>
    friend std::optional<thermal::heat_equation_solution_2d<T>> thermal_stationary_2d(
        std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem);

    template<std::floating_point T, std::signed_integral I>
    friend std::optional<mechanical::mechanical_solution_2d<T>> mechanical_2d(
        std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem, const std::vector<T>& delta_temperature);

public:
    template<std::floating_point T, std::signed_integral I>
    friend void determine_problem(const nlohmann::json& config);
};

template<std::floating_point T, std::signed_integral I>
void problems_1d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    if (parallel::MPI_rank() != 0) {
        logger::warning() << "Calculations are available only on the master process. "
                             "The current process has completed its work." << std::endl;
        return;
    }
    config::check_required_fields(config, {"boundaries", "materials"});
    config::check_optional_fields(config, {"mesh", "auxiliary"});
    if (task.problem == config::problem_t::THERMAL)
        thermal::solve_thermal_1d_problem<T, I>(config, save, task.time_dependency);
    else throw std::domain_error{"Unsupported task. In the one-dimensional case, the following problems are available: \"thermal\""};
}

template<std::floating_point T, std::signed_integral I>
std::optional<thermal::heat_equation_solution_2d<T>> thermal_stationary_2d(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem) {
    using DP = _determine_problem;
    if (!DP::is_thermal(problem))
        return std::nullopt;
    const auto materials = config::thermal_materials_2d<T>{config["materials"], "materials"};
    mesh->neighbours(find_neighbours(*mesh, get_search_radii(materials)));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::NO, !DP::ONLY_LOCAL, DP::SYMMTERIC);
    const auto boundaries_field = problem == config::problem_t::THERMAL ? "boundaries" : "thermal_boundaries";
    return thermal::solve_thermal_2d_problem<T, I>(mesh, materials,
        config::read_thermal_boundary_conditions_2d<T>(config[boundaries_field], boundaries_field),
        config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"}
    );
}

template<std::floating_point T, std::signed_integral I>
void thermal_nonstationary_2d(std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config,
                              const config::save_data& save, const config::problem_t problem) {
    using DP = _determine_problem;
    if (problem != config::problem_t::THERMAL)
        throw std::domain_error{"Mechanical problem does not support time dependence."};
    const auto materials = config::thermal_materials_2d<T>{config["materials"], "materials"};
    mesh->neighbours(find_neighbours(*mesh, get_search_radii(materials)));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::MEMORY, !DP::ONLY_LOCAL, DP::SYMMTERIC);
    thermal::solve_thermal_2d_problem<T, I>(mesh, materials,
        config::read_thermal_boundary_conditions_2d<T>(config["boundaries"], "boundaries"),
        config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"},
        config::time_data<T>{config["time"], "time"},
        save);
}

template<std::floating_point T, std::signed_integral I>
std::optional<mechanical::mechanical_solution_2d<T>> mechanical_2d(
    std::shared_ptr<mesh::mesh_2d<T>>& mesh, const nlohmann::json& config, const config::problem_t problem, const std::vector<T>& delta_temperature) {
    using DP = _determine_problem;
    if (!DP::is_mechanical(problem))
        return std::nullopt;
    const auto materials = config::mechanical_materials_2d<T>{config["materials"], "materials"};
    mesh->neighbours(find_neighbours(*mesh, get_search_radii(materials)));
    mesh::utils::balancing(*mesh, mesh::utils::balancing_t::MEMORY, !DP::ONLY_LOCAL, DP::SYMMTERIC);
    const auto boundaries_field = problem == config::problem_t::MECHANICAL ? "boundaries" : "mechanical_boundaries";
    return mechanical::solve_mechanical_2d_problem<T, I>(mesh, materials,
        config::mechanical_boundaries_conditions_2d<T>{config[boundaries_field], boundaries_field},
        delta_temperature
    );
}

template<std::floating_point T, std::signed_integral I>
void problems_2d(const nlohmann::json& config, const config::save_data& save, const config::task_data& task) {
    if (task.problem == config::problem_t::UNKNOWN)
        throw std::domain_error{"Unknown task. In the two-dimensional case, the following problems are available: "
                                "\"thermal\", \"mechanical\" and \"thermomechanical\""};
    
    using DP = _determine_problem;
    config::check_required_fields(config, DP::get_required_fields(task));
    config::check_optional_fields(config, {"auxiliary"});
    auto mesh = std::make_shared<mesh::mesh_2d<T>>(config::mesh_data<2>{config["mesh"], "mesh"}.path);
    if (task.time_dependency)
        thermal_nonstationary_2d<T, I>(mesh, config, save, task.problem);
    else {
        const std::optional<thermal::heat_equation_solution_2d<T>> thermal_solution = thermal_stationary_2d<T, I>(mesh, config, task.problem);
        const std::optional<mechanical::mechanical_solution_2d<T>> mechanical_solution =
            mechanical_2d<T, I>(mesh, config, task.problem, thermal_solution ? thermal_solution->temperature() : std::vector<T>{});
        save_csv(thermal_solution, mechanical_solution, save);
        save_vtk(thermal_solution, mechanical_solution, save);
    }
}

template<std::floating_point T, std::signed_integral I>
void determine_problem(const nlohmann::json& config) {
    config::check_required_fields(config, {"task"});

    const bool contains_save = config.contains("save");
    const auto save = contains_save ? config::save_data{config["save"], "save"} : config::save_data{};
    if (contains_save)
        _determine_problem::init_save_data(save, config);
    else
        logger::warning() << "There is no \"save\" field in the config. Required data may not be saved." << std::endl;

    if (const config::task_data task{config["task"], "task"}; task.dimension == 1)
        problems_1d<T, I>(config, save, task);
    else if (task.dimension == 2)
        problems_2d<T, I>(config, save, task);
    else throw std::domain_error{"Problem dimension " + std::to_string(task.dimension) + " is not supported"};
}

}
#ifndef NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP
#define NONLOCFEM_ONE_DIMENSIONAM_PROBLEMS_HPP

#include "thermal_problems_utils.hpp"

namespace nonlocal {

template<class T, class I>
void one_dimensional_problems(const nlohmann::json& config, const config::save_data& save, const config::problem_t problem) {
    static const std::set<config::problem_t> available_problems = {
        config::problem_t::THERMAL_STATIONARY,
        config::problem_t::THERMAL_NONSTATIONARY
    };
    if (!available_problems.contains(problem)) {
        throw std::domain_error{
            "In the one-dimensional case, the following problems are available: " +
            init_available_problems_list(available_problems)
        };
    }

    config::check_required_fields(config, {"materials", "boundaries"});
    config::check_required_fields(config, {"mesh", "auxiliary"});
    if (is_thermal_problem(problem)) {
        const config::thermal_materials_1d<T> materials{config["materials"], "materials"};
        const auto mesh = make_mesh_1d(get_segments_data(materials), 
                                       config::mesh_data<1u>{config.value("mesh", nlohmann::json::object()), "mesh"});
        const auto parameters = thermal::make_thermal_parameters<T>(materials.materials);
        const auto auxiliary = config::thermal_auxiliary_data<T>{config.value("auxiliary", nlohmann::json::object()), "auxiliary"};
        const auto boundaries_conditions = thermal::make_thermal_boundaries_conditions_1d(
            config::thermal_boundaries_conditions_1d<T>{config["boundaries"], "boundaries"}
        );

        if (problem == config::problem_t::THERMAL_STATIONARY) {
            auto solution = thermal::stationary_heat_equation_solver_1d<T, I>(
                mesh, parameters, boundaries_conditions,
                thermal::stationary_equation_parameters_1d<T>{
                    .right_part = [value = auxiliary.right_part](const T x) constexpr noexcept { return value; },
                    .initial_distribution = [value = auxiliary.initial_distribution](const T x) constexpr noexcept { return value; },
                    .energy = auxiliary.energy
                }
            );
            if (save.contains("temperature"))
                mesh::utils::save_as_csv(*mesh, solution.temperature(), save.path("temperature", "csv"), save.precision());
            if (save.contains("flux"))
                mesh::utils::save_as_csv(*mesh, solution.calc_flux(), save.path("flux", "csv"), save.precision());
        } else if (problem == config::problem_t::THERMAL_NONSTATIONARY)
            std::cout << "thermal_nonstationary_1d" << std::endl;
    }
    
}
    
}

#endif
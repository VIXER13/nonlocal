#include "make_mesh.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"

namespace {

using T = double;
using I = int64_t;

template<class T>
void save_step(nonlocal::thermal::heat_equation_solution_1d<T>&& solution, const nonlocal::config::save_data& save, const uint64_t step) {
    if (!std::filesystem::exists(save.folder()))
        std::filesystem::create_directories(save.folder());
    if (save.contains("temperature"))
        nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.temperature(), save.path("", ".csv", std::to_string(step) + "temperature"));
    if (save.contains("flux"))
        nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.calc_flux(), save.path("", ".csv", std::to_string(step) + "flux"));
}

}

int main(const int argc, const char *const *const argv) {
    if (argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        const nonlocal::config::nonstationary_thermal_1d_data<T> config_data{
            nonlocal::config::read_json(std::filesystem::path{argv[1]})
        };

        const auto mesh = nonlocal::make_mesh(config_data.materials, config_data.element_order, config_data.quadrature_order);
        const auto parameters = nonlocal::make_thermal_parameters(config_data.materials);
        const nonlocal::thermal::thermal_boundaries_conditions_1d<T> boundaries_conditions =  {
            nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("left")),
            nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("right"))
        };

        nonlocal::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, config_data.nonstationary.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = config_data.equation.initial_distribution](const T x) constexpr noexcept { return init_dist; });
        save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, config_data.save, 0u);
        for(const uint64_t step : std::ranges::iota_view{1u, config_data.nonstationary.steps_cont + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = config_data.equation.right_part](const T x) constexpr noexcept { return right_part; });
            if (step % config_data.nonstationary.save_frequency == 0) {
                std::cout << "step = " << step << std::endl;
                save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, config_data.save, step);
            }
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
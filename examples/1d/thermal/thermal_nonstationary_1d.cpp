#include "make_mesh.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"

namespace {

template<class T>
void save_step(nonlocal::thermal::heat_equation_solution_1d<T>&& solution, const nonlocal::config::save_data& save, const uint64_t step) {
    std::cout << "step = " << step << std::endl;
    std::cout << "integral = " << nonlocal::mesh::utils::integrate(solution.mesh(), solution.temperature()) << std::endl;
    const auto save_vector = [&solution, &save, step](const std::vector<T>& x, const std::string& name) {
        const std::filesystem::path path = save.make_path(std::to_string(step) + save.get_name(name), ".csv");
        nonlocal::mesh::utils::save_as_csv(solution.mesh(), x, path, save.precision());
    };
    if (save.contains("temperature"))
        save_vector(solution.temperature(), "temperature");
    if (save.contains("flux"))
        save_vector(solution.calc_flux(), "flux");
}

}

int main(const int argc, const char *const *const argv) {
    if (argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        using T = double;
        using I = int64_t;
        const nonlocal::config::nonstationary_thermal_data<T, 1> config_data{nonlocal::config::read_json(std::filesystem::path{argv[1]})};
        std::cout.precision(config_data.other.get("precision", std::cout.precision()).asInt());

        const auto mesh = nonlocal::make_mesh(config_data.materials, config_data.elements.element_order, config_data.elements.quadrature_order);
        const auto parameters = nonlocal::make_thermal_parameters(config_data.materials);
        const nonlocal::thermal::thermal_boundaries_conditions_1d<T> boundaries_conditions =  {
            nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("left")),
            nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("right"))
        };
        nonlocal::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, config_data.nonstationary.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = config_data.equation.initial_distribution](const T x) constexpr noexcept { return init_dist; });
            
        if (!std::filesystem::exists(config_data.save.folder()))
            std::filesystem::create_directories(config_data.save.folder());
        if (config_data.save.contains("config"))
            nonlocal::config::save_json(config_data.save.path("config", ".json"), config_data.to_json());
        save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, config_data.save, 0u);
        for(const uint64_t step : std::ranges::iota_view{1u, config_data.nonstationary.steps_cont + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = config_data.equation.right_part](const T x) constexpr noexcept { return right_part; });
            if (step % config_data.nonstationary.save_frequency == 0)
                save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, config_data.save, step);
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
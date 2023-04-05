#include "make_mesh.hpp"
#include "thermal/nonstationary_heat_equation_solver_2d.hpp"

namespace {

template<class T, class I>
void save_step(nonlocal::thermal::heat_equation_solution_2d<T, I>&& solution, const nonlocal::config::save_data& save, const uint64_t step) {
    if (parallel_utils::MPI_rank() != 0)
        return;

    std::cout << "step = " << step << std::endl;
    const auto save_vector = [&solution, &save, step](const std::vector<T>& x, const std::string& name) {
        const std::filesystem::path path = save.make_path(std::to_string(step) + save.get_name(name), ".csv");
        nonlocal::mesh::utils::save_as_csv(path, solution.mesh().container(), x, save.precision());
    };

    solution.calc_flux();
    if (save.contains("temperature"))
        save_vector(solution.temperature(), "temperature");
    if (save.contains("flux_x"))
        save_vector(solution.flux()[0], "flux_x");
    if (save.contains("flux_y"))
        save_vector(solution.flux()[1], "flux_y");
    if (save.contains("vtk"))
        solution.save_as_vtk(save.make_path(std::to_string(step) + save.get_name("vtk"), ".vtk"));
}

}

int main(const int argc, const char *const *const argv) {
    if(argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

#ifdef MPI_USED
    MPI_Init(&argc, &argv);
#endif

    int exit_code = EXIT_SUCCESS;
    try {
        using T = double;
        using I = int64_t;
        const nonlocal::config::nonstationary_thermal_data<T, 2> config_data{nonlocal::config::read_json(std::filesystem::path{argv[1]})};
        std::cout.precision(config_data.other.get("precision", std::cout.precision()).asInt());

        const auto mesh = nonlocal::make_mesh<T, I>(config_data.mesh.mesh, config_data.materials);
        const auto parameters = nonlocal::make_parameters<T>(config_data.materials);
        const auto boundaries_conditions = nonlocal::make_boundaries_conditions(config_data.boundaries);
        nonlocal::thermal::nonstationary_heat_equation_solver_2d<T, I, I> solver{mesh, config_data.nonstationary.time_step};
        solver.compute(parameters, boundaries_conditions,
            [init_dist = config_data.equation.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return init_dist; });

        if (!std::filesystem::exists(config_data.save.folder()))
            std::filesystem::create_directories(config_data.save.folder());
        if (config_data.save.contains("config"))
            nonlocal::config::save_json(config_data.save.path("config", ".json"), config_data.to_json());
        save_step(nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, parameters, solver.temperature()}, config_data.save, 0u);
        for(const uint64_t step : std::ranges::iota_view{1u, config_data.nonstationary.steps_cont + 1}) {
            solver.calc_step(boundaries_conditions,
                [right_part = config_data.equation.right_part](const std::array<T, 2>& x) constexpr noexcept { return right_part; });
            if (step % config_data.nonstationary.save_frequency == 0)
                save_step(nonlocal::thermal::heat_equation_solution_2d<T, I>{mesh, parameters, solver.temperature()}, config_data.save, step);
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit_code = EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        exit_code = EXIT_FAILURE;
    }

#ifdef MPI_USED
    MPI_Finalize();
#endif

    return exit_code;
}
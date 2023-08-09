#include "make_mesh.hpp"

namespace {

template<class T, class I>
void save_data(const nonlocal::thermal::heat_equation_solution_2d<T, I>& solution, 
               const nonlocal::config::stationary_thermal_data<T, 2>& config_data) {
    const auto& [TX, TY] = solution.flux();
    if (!std::filesystem::exists(config_data.save.folder()))
        std::filesystem::create_directories(config_data.save.folder());
    if (config_data.save.contains("config"))
        nonlocal::config::save_json(config_data.save.path("config", ".json"), config_data);
    if (config_data.save.contains("temperature"))
        nonlocal::mesh::utils::save_as_csv(config_data.save.path("temperature", ".csv"), solution.mesh().container(), solution.temperature(), config_data.save.precision());
    if (config_data.save.contains("flux_x"))
        nonlocal::mesh::utils::save_as_csv(config_data.save.path("flux_x", ".csv"), solution.mesh().container(), TX, config_data.save.precision());
    if (config_data.save.contains("flux_y"))
        nonlocal::mesh::utils::save_as_csv(config_data.save.path("flux_y", ".csv"), solution.mesh().container(), TY, config_data.save.precision());
    if (config_data.save.contains("vtk"))
        solution.save_as_vtk(config_data.save.path("vtk", ".vtk"));
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
        const nonlocal::config::stationary_thermal_data<T, 2> config_data{nonlocal::config::read_json(std::filesystem::path{argv[1]})};
        std::cout.precision(config_data.other.get("precision", std::cout.precision()).asInt());

        constexpr T sigma = T{7};
        constexpr T epsilon = T{10.0};

        auto parameters = nonlocal::make_parameters<T>(config_data.materials);
        parameters["DEFAULT"].physical = std::make_shared<nonlocal::thermal::parameter_2d<T, nonlocal::coefficients_t::SOLUTION_DEPENDENT>>(
            [epsilon, sigma](const std::array<T, 2>& x, const T solution) { return (T{0.1} + epsilon * std::pow(std::abs(solution), sigma)); }
            //[](const std::array<T, 2>& x) { return 1.0; }
        );
        const auto additional_parameters = nonlocal::thermal::stationary_equation_parameters_2d<T>{
            .right_part = [value = config_data.equation.right_part](const std::array<T, 2>& x) constexpr noexcept { return value; },
            .initial_distribution = [value = config_data.equation.initial_distribution](const std::array<T, 2>& x) constexpr noexcept { return value; },
            .energy = config_data.equation.energy
        };

        const auto mesh = nonlocal::make_mesh<T, I>(config_data.mesh.path, config_data.materials);
        auto solution = nonlocal::thermal::stationary_heat_equation_solver_nonlinear_2d<I>(
            mesh, parameters, 
            nonlocal::make_boundaries_conditions(config_data.boundaries), 
            additional_parameters
        );
        solution.calc_flux();

        if (parallel_utils::MPI_rank() == 0)
            save_data(solution, config_data);
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
#include "make_mesh.hpp"

int main(const int argc, const char *const *const argv) {
    if (argc != 2) {
        std::cerr << "Input format [program name] <path/to/config.json>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        using T = double;
        using I = int64_t;
        const nonlocal::config::stationary_thermal_data<T, 1> config_data{nonlocal::config::read_json(std::filesystem::path{argv[1]})};
        std::cout.precision(config_data.other.get("precision", std::cout.precision()).asInt());

        constexpr T eps = T{10};
        constexpr T sigma = T{1};

        auto params = nonlocal::make_thermal_parameters(config_data.materials);
        params[0].physical = std::make_shared<nonlocal::thermal::parameter_1d<T, nonlocal::coefficients_t::SOLUTION_DEPENDENT>>(
                // [eps] (const T x, const T temp) { return (T{1} + eps * temp); },
                [sigma] (const T x, const T temp) { return std::pow(temp, sigma); },
                //[] (const T x, const T temp) { return lambda_0 * (1 + eps * temp); },
                config_data.materials[0].physical.capacity,
                config_data.materials[0].physical.density
            );


        const auto mesh = nonlocal::make_mesh(config_data.materials, config_data.mesh.element_order, config_data.mesh.quadrature_order);


        // auto solution = nonlocal::thermal::stationary_heat_equation_solver_1d<T, I>(
        //     mesh, params, //nonlocal::make_thermal_parameters(config_data.materials),
        //     nonlocal::thermal::thermal_boundaries_conditions_1d<T>{
        //         nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("left")),
        //         nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("right"))
        //     },
        //     [value = config_data.equation.right_part](const T x) constexpr noexcept { return value; },
        //     config_data.equation.energy
        // );

        auto solution = nonlocal::thermal::stationary_nonlinear_heat_equation_solver_1d<T, I>(
            mesh, params, //nonlocal::make_thermal_parameters(config_data.materials),
            nonlocal::thermal::thermal_boundaries_conditions_1d<T>{
                nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("left")),
                nonlocal::make_boundary_condition<T>(config_data.boundaries.conditions.at("right"))
            },
            [value = config_data.equation.right_part](const T x) constexpr noexcept { return value; },
            config_data.equation.initial_distribution,
            config_data.equation.energy
        );

        
        std::cout << "integral = " << nonlocal::mesh::utils::integrate(*mesh, solution.temperature()) << std::endl;
        if (!std::filesystem::exists(config_data.save.folder()))
            std::filesystem::create_directories(config_data.save.folder());
        if (config_data.save.contains("config"))
            nonlocal::config::save_json(config_data.save.path("config", ".json"), config_data);
        if (config_data.save.contains("temperature"))
            nonlocal::mesh::utils::save_as_csv(*mesh, solution.temperature(), config_data.save.path("temperature", ".csv"), config_data.save.precision());
        if (config_data.save.contains("flux"))
            nonlocal::mesh::utils::save_as_csv(*mesh, solution.calc_flux(), config_data.save.path("flux", ".csv"), config_data.save.precision());
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
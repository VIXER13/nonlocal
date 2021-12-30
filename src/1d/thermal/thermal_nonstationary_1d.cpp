#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver.hpp"
#include "influence_functions_1d.hpp"
#include <iostream>

int main(const int argc, const char *const *const argv) {
    if (argc < 6) {
        std::cerr << "run format: program_name <element_type> <elements_count> <p1> <r> <save_path>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<double>>(
            nonlocal::make_element<double>(nonlocal::element_type(std::stoi(argv[1]))),
            std::stoull(argv[2]), std::array{0., 1.});
        const nonlocal::nonlocal_parameters_1d<double> nonloc_parameters = {
            .p1 = std::stod(argv[3]),
            .r  = std::stod(argv[4])
        };
        const nonlocal::thermal::heat_equation_parameters_1d<double> equation_parameters = {
            .lambda = 1,
            .alpha = {2., 2. / 19.}
        };
        const nonlocal::nonstationary_solver_parameters_1d<double> sol_parameters = {
            .save_path = argv[5],
            .time_interval = {0, 1},
            .steps = 100,
            .save_freq = 1
        };

        nonlocal::thermal::nonstationary_heat_equation_solver_1d<double, int>(
            nonloc_parameters, equation_parameters, sol_parameters, mesh,
            {
                nonlocal::thermal::boundary_condition_t::FLUX, [](const double) constexpr noexcept { return  1; },
                nonlocal::thermal::boundary_condition_t::FLUX, [](const double) constexpr noexcept { return -1; }
            },
            [](const double) { return 0; },
            [](const double t, const double x) { return 0; },
            nonlocal::influence::polynomial_1d<double, 2, 1>{nonloc_parameters.r}
        );
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
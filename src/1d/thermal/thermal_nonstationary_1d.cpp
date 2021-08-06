#include <iostream>
#include <ostream>
#include "heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "make_element.hpp"

namespace {

constexpr uintmax_t factorial(uintmax_t n) {
    return n > 0 ? n * factorial(n - 1) : 1;
}

template<class T, uintmax_t M>
T impulse(const T t) noexcept {
    using metamath::function::power;
    static constexpr T coeff = power<M>(T{M}) / factorial(M);
    return coeff * power<M>(t) * std::exp(-t * M);
}

}

int main(int argc, char** argv) {
    if (argc < 8) {
        std::cerr << "run format: program_name <element_type> <elements_count> <p1> <r> <save_path> <steps> <save_frequence>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<double>>(
            nonlocal::make_element<double>(nonlocal::element_type(std::stoi(argv[1]))),
            std::stoull(argv[2]), std::array{0., 5.});

        nonlocal::heat::equation_parameters<double> parameters;
        parameters.p1 = std::stod(argv[3]);
        parameters.r = std::stod(argv[4]);
        mesh->calc_neighbours_count(parameters.r);
        nonlocal::heat::heat_equation_solver_1d<double> solver{mesh};

        nonlocal::heat::solver_parameters<double> sol_parameters;
        sol_parameters.save_path = argv[5];
        sol_parameters.time_interval[0] = 0;
        sol_parameters.time_interval[1] = 10;
        sol_parameters.steps = std::stoull(argv[6]);
        sol_parameters.save_freq = std::stoull(argv[7]);

        solver.nonstationary(sol_parameters, parameters,
             {
                 std::pair{
                     nonlocal::heat::FLOW,
                     [](const double t) noexcept { return impulse<double, 2>(t); }
                 },
                 std::pair{
                     nonlocal::heat::FLOW,
                     [](const double t) noexcept { return 0; }
                 }
             },
             [](const double x) noexcept { return 0; },
             [](const double t, const double x) noexcept { return 0; },
             nonlocal::influence::polynomial_1d<double, 2, 1>{parameters.r}
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
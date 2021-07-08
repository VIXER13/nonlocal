#include <iostream>
#include <ostream>
#include "heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "make_element.hpp"

namespace {

template<class Vector>
void save_as_csv(const std::string& path, const Vector& x) {
    std::ofstream csv{path};
    //csv.precision(std::numeric_limits<T>::max_digits10);
    const double h = 1. / (x.size() - 1);
    for(size_t i = 0; i < x.size(); ++i)
        csv << i * h << ',' << x[i] << '\n';
}

}

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "run format: program_name <element_type> <elements_count> <section>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<mesh::mesh_1d<double>>(
            nonlocal::make_element<double>(nonlocal::element_type(std::stoi(argv[1]))),
            std::stoull(argv[2]), std::array{std::stod(argv[3]), std::stod(argv[4])}
        );

        nonlocal::solver_parameters<double> sol_parameters;
        sol_parameters.save_path = "./impulseNonLocP12R05/";
        //sol_parameters.save_path = "./impulseLoc/";
        sol_parameters.time_interval[0] = 0;
        sol_parameters.time_interval[1] = 10;
        sol_parameters.steps = 10000;
        sol_parameters.save_freq = 1;

        nonlocal::heat::equation_parameters<double> parameters;
        parameters.p1 = 0.5;
        parameters.r = 0.5;
        mesh->calc_neighbours_count(parameters.r);
        nonlocal::heat::heat_equation_solver_1d<double> solver{mesh};
        solver.nonstationary(sol_parameters, parameters,
             {
                 std::pair{
                     nonlocal::boundary_condition_t::SECOND_KIND,
                     [](const double t) { return 4 * t * t * std::exp(-2 * t); }
                 },
                 std::pair{
                     nonlocal::boundary_condition_t::SECOND_KIND,
                     [](const double t) { return 0; }
                 }
             },
             [](const double x) { return 0; },
             [](const double t, const double x) { return 0; },
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
#include <iostream>
#include <ostream>
#include "heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "make_element.hpp"

namespace {

template<class T, class Vector>
void save_as_csv(const std::string& path, const Vector& x, const std::array<T, 2>& section) {
    std::ofstream csv{path};
    const T h = (section.back() - section.front()) / (x.size() - 1);
    for(size_t i = 0; i < x.size(); ++i)
        csv << section.front() +  i * h << ',' << x[i] << '\n';
}

}

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "run format: program_name <element_type> <elements_count> <p1> <r> <save_name>" << std::endl;
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<mesh::mesh_1d<double>>(
            nonlocal::make_element<double>(nonlocal::element_type(std::stoi(argv[1]))),
            std::stoull(argv[2]), std::array{0., 1.});

        nonlocal::heat::equation_parameters<double> parameters;
        parameters.p1 = std::stod(argv[3]);
        parameters.r = std::stod(argv[4]);
        mesh->calc_neighbours_count(parameters.r);
        nonlocal::heat::heat_equation_solver_1d<double> solver{mesh};

        auto solution = solver.stationary(
            parameters,
            {
                std::pair{nonlocal::boundary_condition_t::SECOND_KIND, -1.},
                std::pair{nonlocal::boundary_condition_t::SECOND_KIND,  1.},
            },
            [](const double x) { return 0; },
            nonlocal::influence::polynomial_1d<double, 2, 1>{parameters.r}
        );
        save_as_csv(argv[5], solution, mesh->section());
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
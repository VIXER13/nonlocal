#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include <iostream>
#include <fstream>

namespace {

template<class Vector>
void save_step(const std::string& path, const Vector& vec, const nonlocal::mesh::mesh_1d<double>& mesh, const uintmax_t step) {
    using namespace std::literals;
    std::ofstream csv{path + "/"s + std::to_string(step) + ".csv"};
    csv.precision(std::numeric_limits<double>::max_digits10);
    const double h = (mesh.section().back() - mesh.section().front()) / (mesh.nodes_count() - 1);
    for(const size_t i : std::views::iota(size_t{0}, mesh.nodes_count()))
        csv << mesh.section().front() + i * h << ',' << vec[i] << '\n';
}

}

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
        const double p1 = std::stod(argv[3]),
                     r  = std::stod(argv[4]);
        const double tau = 0.01;
        const nonlocal::thermal::heat_equation_parameters_1d<double> equation_parameters = {
            .lambda = 1,
            .alpha = {2., 2. / 19.}
        };

        const std::array<nonlocal::nonstatinary_boundary_1d_t<nonlocal::thermal::boundary_condition_t, double>, 2>
            boundary_condition = {
                nonlocal::thermal::boundary_condition_t::FLUX, [](const double) constexpr noexcept { return  1; },
                nonlocal::thermal::boundary_condition_t::FLUX, [](const double) constexpr noexcept { return -1; }
        };

        mesh->calc_neighbours_count(r);
        nonlocal::thermal::nonstationary_heat_equation_solver_1d<double, int> solver{mesh, tau};
        solver.compute(equation_parameters,
            nonlocal::boundary_type(boundary_condition),
            [](const double) constexpr noexcept { return 0; },
            p1, nonlocal::influence::polynomial_1d<double, 2, 1>{r}
        );
        save_step(argv[5], solver.temperature(), *mesh, 0);
        for(const uintmax_t step : std::views::iota(1, 101)) {
            solver.calc_step(boundary_condition, [](const double t, const double x) constexpr noexcept { return 0; });
            save_step(argv[5], solver.temperature(), *mesh, step);
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
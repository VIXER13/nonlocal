#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include <iostream>

namespace {

using T = double;
using I = int64_t;

void save_step(nonlocal::thermal::heat_equation_solution_1d<T>&& solution, const std::filesystem::path& folder, const uintmax_t step) {
    using namespace std::literals;
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.temperature(), folder / ("T" + std::to_string(step) + ".csv"));
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.calc_flux(), folder / ("Flux" + std::to_string(step) + ".csv"));
}

}

int main(const int argc, const char *const *const argv) {
    try {
        std::cout.precision(3);
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
            nonlocal::make_element<T>(nonlocal::element_type::QUADRATIC),
            std::vector<nonlocal::mesh::segment_data<T>>{
                {.length = 0.15, .elements = 100},
                {.length = 0.25, .elements = 100},
                {.length = 0.35, .elements = 100},
                {.length = 0.15, .elements = 100}
            }
        );
        const std::vector<T> radii = {
            0.05, 
            0., 
            0.1, 
            0.
        };
        const T p1 = 0.5;
        const T tau = 0.001;
        std::vector<nonlocal::equation_parameters<1, T, nonlocal::thermal::parameters_1d>> parameters;
        for(const auto [conductivity, radius, local_weight] : {std::tuple{ 1., radii[0], p1 }, 
                                                               std::tuple{ 7., radii[1], 1. },
                                                               std::tuple{ 3., radii[2], p1 },
                                                               std::tuple{10., radii[3], 1. }
                                                               }) {
            parameters.push_back({
                .physical = {
                    .conductivity = conductivity
                },
                .model = {
                    .influence = nonlocal::influence::polynomial_1d<T, 1, 1>{radius},
                    .local_weight = local_weight
                }
            });
        }
        if (nonlocal::theory_type(p1) == nonlocal::theory_t::NONLOCAL)
            mesh->find_neighbours(radii);

        const std::array<std::unique_ptr<nonlocal::thermal::stationary_thermal_boundary_condition_1d>, 2> boundary_condition = {
            std::make_unique<nonlocal::thermal::stationary_flux_1d<T>>(1.),
            std::make_unique<nonlocal::thermal::stationary_flux_1d<T>>(-1.)
        };

        nonlocal::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, tau};
        solver.compute(
            parameters,
            boundary_condition,
            [](const double) constexpr noexcept { return 0; }
        );
        save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, "./sol", 0);
        for(const uintmax_t step : std::ranges::iota_view{1u, 101u}) {
            solver.calc_step(boundary_condition, [](const double x) constexpr noexcept { return 0; });
            save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, "./sol", step);
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
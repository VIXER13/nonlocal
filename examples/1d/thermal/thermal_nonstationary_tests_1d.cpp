#include "make_element.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"
#include "metamath.hpp"

#include <iostream>
#include <cmath>

namespace {

using T = double;
using I = int64_t;

void save_step(nonlocal::thermal::heat_equation_solution_1d<T>&& solution, const std::filesystem::path& folder, const uintmax_t step) {
    using namespace std::literals;
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.temperature(), folder / ("T" + std::to_string(step) + ".csv"));
    nonlocal::mesh::utils::save_as_csv(solution.mesh(), solution.calc_flux(), folder / ("Flux" + std::to_string(step) + ".csv"));
}

void save_info(const std::vector<nonlocal::equation_parameters<1, T, nonlocal::thermal::parameters_1d>>& parameters,
               const std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2>& boundary_condition,
               const T tau, const std::vector<T>& radii, const std::vector<nonlocal::mesh::segment_data<T>>& segment_data,
               const std::filesystem::path& folder) {
    std::ofstream info_file;
    const std::string way_to_file = folder /  ("calculation_info.txt");
    std::cout << "Info about calculation will be writen in " << way_to_file << "\n";
    info_file.open(way_to_file);

    info_file << "tau = "                << tau               << "\n";
    info_file << "number of segments = " << parameters.size() << "\n";
    for (size_t i = 0; i < parameters.size(); ++i) {
    info_file << "---------segment-№"<< i + 1 << "-----------------\n";
    info_file << "length = "             << std::to_string(segment_data[i].length)              << "\n";
    info_file << "number of elements = " << std::to_string(segment_data[i].elements)            << "\n";
    info_file << "conductivity = "       << std::to_string(parameters[i].physical.conductivity) << "\n";
    info_file << "density = "            << std::to_string(parameters[i].physical.density)      << "\n";
    info_file << "capacity = "           << std::to_string(parameters[i].physical.capacity)     << "\n";
    info_file << "local_weight (p1) = "  << std::to_string(parameters[i].model.local_weight)    << "\n";
    info_file << "radii = "              << std::to_string(radii[i])                            << "\n";
    info_file << "-------------------------------------\n";
    }

    info_file.close();
}

}

int main(const int argc, const char *const *const argv) {
    try {
        std::cout.precision(3);
        const std::vector<nonlocal::mesh::segment_data<T>> segment_data = {
                {.length = 1., .elements = 100}
                //{.length = 0.25, .elements = 100},
                //{.length = 0.35, .elements = 100},
                //{.length = 0.15, .elements = 100}
            };
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
            nonlocal::make_element<T>(nonlocal::element_type::LINEAR),
            segment_data 
        );
        const std::vector<T> radii = {
            0.05
            //0., 
            //0.1, 
            //0.
        };
        //const T p1 = 0.5;
        const T p1 = 1.;
        const T tau = 0.01;
        std::vector<nonlocal::equation_parameters<1, T, nonlocal::thermal::parameters_1d>> parameters;
        for(const auto [conductivity, radius, local_weight] : {std::tuple{ 1., radii[0], p1 }
                                                            // std::tuple{ 7., radii[1], 1. },
                                                            // std::tuple{ 3., radii[2], p1 },
                                                            // std::tuple{10., radii[3], 1. }
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

        const size_t segment = 0;

        //Initial temperature portrait
        auto initial_dist = [](const double ) constexpr noexcept { return 0.; };

        //Flux BC right
        auto flux = [](const double t ) {
         return 0; 
        };

        //Convection BC left
        const double heat_transfer = 10.;
        auto ambient_temperature = [](const double t ) { 
            return 10; 
        };



        nonlocal::thermal::nonstationary_heat_equation_solver_1d<T, I> solver{mesh, tau};
        const std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2> boundary_condition = {
            std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(0, heat_transfer, ambient_temperature(0), 0.0, 0),
            std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(0, heat_transfer, ambient_temperature(0), 0.0, 0)
        };


        solver.compute(
            parameters,
            boundary_condition,
            initial_dist
        );

        save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, argv[1], 0);
        //save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, "./sol", 0);
        for(const uintmax_t step : std::ranges::iota_view{1u, 1001u}) {
            const T time = step * tau;

            const std::array<std::unique_ptr<nonlocal::thermal::thermal_boundary_condition_1d<T>>, 2> boundary_condition = {
                std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(0, heat_transfer, ambient_temperature(0), 0.0, solver.temperature()[0]),
                std::make_unique<nonlocal::thermal::combined_flux_1d<T>>(0, heat_transfer, ambient_temperature(0), 0.0, solver.temperature()[solver.temperature().size() - 1])
            };

            //Right hand side
            auto rhs = [](const double x) { 
                return 0; 
            };

            solver.calc_step(boundary_condition, rhs);
            save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, argv[1], step);
            //save_step(nonlocal::thermal::heat_equation_solution_1d<T>{mesh, parameters, solver.temperature()}, "./sol", step);
        }
        save_info(parameters, boundary_condition, tau, radii, segment_data, argv[1]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
#include "make_element.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
//#include "influence_functions_1d.hpp"
#include "mesh_1d_utils.hpp"
#include <iostream>
#include <fstream>

namespace {

using T = double;
using I = int64_t;

template<class Vector>
void save_as_csv(const std::string& path, const Vector& x, const std::array<T, 2>& section) {
    std::ofstream csv{path};
    const T h = (section.back() - section.front()) / (x.size() - 1);
    for(size_t i = 0; i < x.size(); ++i)
        csv << section.front() +  i * h << ',' << x[i] << '\n';
}

}

int main(const int argc, const char *const *const argv) {
    //if (argc < 6) {
    //    std::cerr << "run format: program_name <element_type> <elements_count> <p1> <r> <save_name> {grad_name}" << std::endl;
    //    return EXIT_FAILURE;
    //}

    try {
        std::cout.precision(3);
        const auto mesh = std::make_shared<nonlocal::mesh::mesh_1d<T>>(
            nonlocal::make_element<T>(nonlocal::element_type::QUADRATIC),
            std::vector{
                nonlocal::mesh::segment_data{.length = 0.25, .elements = 100},
                nonlocal::mesh::segment_data{.length = 0.25, .elements = 100},
                nonlocal::mesh::segment_data{.length = 0.25, .elements = 100},
                nonlocal::mesh::segment_data{.length = 0.25, .elements = 100}
            }
        );

        const std::vector<nonlocal::equation_parameters<1, T, nonlocal::thermal::stationary_equation_parameters_1d>> parameters = {
            {   .physical = {
                    .conductivity = 1
                }
            },
            {   .physical = {
                    .conductivity = 7
                }
            },
            {   .physical = {
                    .conductivity = 3
                }
            },
            {   .physical = {
                    .conductivity = 10
                }
            },
        };

        const auto solution = nonlocal::thermal::stationary_heat_equation_solver_1d<T, I>(
            mesh, parameters,
            {   nonlocal::thermal::boundary_condition_t::TEMPERATURE, 1.,
                nonlocal::thermal::boundary_condition_t::TEMPERATURE, 0.,
            },
            [](const T x) constexpr noexcept { return 0; }
        );
        //std::cout << "integral = " << nonlocal::utils::integrate(*mesh, solution) << std::endl;
        save_as_csv("Tsegmented.csv", solution, {T{0}, T{1}});

        if (argc == 7) {
            //const auto flux = nonlocal::utils::calc_flux(*mesh, solution, p1, influence_function);
            //save_as_csv(argv[6], flux, mesh->section());
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
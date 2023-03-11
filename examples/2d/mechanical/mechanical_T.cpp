#include <iostream>
#include "equilibrium_equation_2d.hpp"
#include "influence_functions_2d.hpp"
#include "mesh_container_2d_utils.hpp"

namespace {

using T = double;
using I = int64_t;

}

int main(const int argc, const char *const *const argv) {
#ifdef MPI_USE
    MPI_Init(&argc, &argv);
#endif

    if(argc < 6) {
        std::cerr << "Input format [program name] <path to mesh> <r1> <r2> <p1> <path_to_save>" << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(10);
        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
        const std::array<T, 2> r = {std::stod(argv[2]), std::stod(argv[3])};
        const T p1 = std::stod(argv[4]);
        nonlocal::mechanical::mechanical_parameters_2d<T> parameters;
        parameters.materials["Material1"] = {
            .model = {
                .influence = nonlocal::influence::polynomial_2d<T, 2, 1>{r},
                .local_weight = T{1}
            },
            .physical = {
                .young_modulus = 1,
                .poissons_ratio = 0.3
            }
        };
        parameters.materials["Material2"] = {
            .model = {
                .influence = nonlocal::influence::normal_distribution_2d<T>{std::stod(argv[2])},
                .local_weight = p1
            },
            .physical = {
                .young_modulus = 2,
                .poissons_ratio = 0.2
            }
        };
        parameters.materials["Material3"] = {
            .model = {
                .influence = nonlocal::influence::polynomial_2d<T, 2, 1>{r},
                .local_weight = T{1}
            },
            .physical = {
                .young_modulus = 3,
                .poissons_ratio = 0.1
            }
        };

        // mesh->find_neighbours({
        //     {"Material2", 3 * std::stod(argv[2])}
        //     //{"Material4", std::stod(argv[2]) + 0.05}
        // });

        nonlocal::mechanical::mechanical_boundaries_conditions_2d<T> boundary_conditions;
        boundary_conditions["Left"] = {
            std::make_unique<nonlocal::mechanical::displacement_2d<T>>(0.),
            std::make_unique<nonlocal::mechanical::displacement_2d<T>>(0.)
        };
        boundary_conditions["Right"] = {
            std::make_unique<nonlocal::mechanical::pressure_2d<T>>(1.),
            std::make_unique<nonlocal::mechanical::pressure_2d<T>>(0.)
        };
        
        static constexpr auto right_part = [](const std::array<T, 2>& x) constexpr noexcept {
            return std::array<T, 2>{};
        };

        auto solution = nonlocal::mechanical::equilibrium_equation<I>(
            mesh, parameters, boundary_conditions, right_part
        );

        if (parallel_utils::MPI_rank() == 0) {
            using namespace std::literals;
            solution.calc_strain_and_stress();
            solution.save_as_vtk(argv[5] + "/mechanical.vtk"s);
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

#ifdef MPI_USE
    MPI_Finalize();
#endif
    return 0;
}
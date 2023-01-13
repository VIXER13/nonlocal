#include <iostream>
#include "stationary_heat_equation_solver_2d.hpp"
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

    if(argc < 5) {
        std::cerr << "Input format [program name] <path to mesh> <r1> <r2> <p1>" << std::endl;
#ifdef MPI_USE
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(3);
        auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
        const std::array<T, 2> r = {std::stod(argv[2]), std::stod(argv[3])};
        const T p1 = std::stod(argv[4]);
        nonlocal::thermal::parameter_2d<T> parameters = {
            .conductivity = {T{1}, T{0},
                             T{0}, T{1}},
            .material = nonlocal::material_t::ISOTROPIC
        };

        if (nonlocal::theory_type(p1) == nonlocal::theory_t::NONLOCAL) {
            mesh->find_neighbours(std::max(r[0], r[1]));
        }

        nonlocal::boundaries_conditions_2d<T, nonlocal::physics_t::THERMAL, 1> boundary_conditions;

        boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(1);
        boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(-1);

        // boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(
        //      T{10}, [](const std::array<T, 2>& x) constexpr noexcept { return std::abs(x[1] - 0.5); }
        // );
        // boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(
        //      T{10}, [](const std::array<T, 2>& x) constexpr noexcept { return std::abs(x[1] - 0.5); }
        // );
        // boundary_conditions["Up"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
        // );
        // boundary_conditions["Down"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        //     [](const std::array<T, 2>& x) constexpr noexcept { return T{0}; }
        // );
        
        static constexpr auto right_part = [](const std::array<T, 2>& x) constexpr noexcept {
            return T{0} ;
        };

        auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
            mesh, parameters, boundary_conditions, right_part, p1, nonlocal::influence::polynomial_2d<T, 1, 1>{r}
        );

        if (parallel_utils::MPI_rank() == 0) {
            const auto& [TX, TY] = solution.calc_flux();
            //std::cout << "Energy = " << solution.calc_energy() << std::endl;
            solution.save_as_vtk("heat.vtk");
            using namespace std::literals;
            nonlocal::mesh::utils::save_as_csv("T.csv", mesh->container(), solution.temperature());
            nonlocal::mesh::utils::save_as_csv("TX.csv", mesh->container(), TX);
            nonlocal::mesh::utils::save_as_csv("TY.csv", mesh->container(), TY);
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
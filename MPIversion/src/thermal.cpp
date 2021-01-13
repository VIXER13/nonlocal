#include <petsc.h>
#include "influence_functions.hpp"
#include "heat_equation_solver.hpp"

namespace {

void save_raw_data(const mesh::mesh_2d<double>& msh, const Eigen::Matrix<double, Eigen::Dynamic, 1>& T) {
    std::ofstream Tout{"T.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i)
        Tout << msh.node(i)[0] << ',' << msh.node(i)[1] << ',' << T[i] << '\n';
}

}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);
        omp_set_num_threads(1);

        static constexpr double r = 0.2, p1 = 1;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};

        nonlocal::heat::heat_equation_solver<double, int> fem_sol{msh};

        const auto T = fem_sol.stationary(
            { // Граничные условия
                {   // Down
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return -1; },
                },

                {   // Right
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                {   // Up
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 1; },
                },

                {   // Left
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                }
            },
            [](const std::array<double, 2>&) { return 0; }, // Правая часть
            r,  // Радиус влияния
            p1, // Вес
            bell // Функция влияния
        );

//        const auto T = fem_sol.stationary(
//            { // Граничные условия
//                {   // Down
//                    nonlocal::heat::boundary_t::FLOW,
//                    [](const std::array<double, 2>& x) { return -1; },
//                },
//
//                {   // Right
//                    nonlocal::heat::boundary_t::FLOW,
//                    [](const std::array<double, 2>& x) { return 0; },
//                },
//
//                {   // Up
//                    nonlocal::heat::boundary_t::FLOW,
//                    [](const std::array<double, 2>& x) { return 1; },
//                },
//
//                {   // Left
//                    nonlocal::heat::boundary_t::FLOW,
//                    [](const std::array<double, 2>& x) { return 0; },
//                }
//            },
//            [](const std::array<double, 2>&) { return -4; }, // Правая часть
//            r,  // Радиус влияния
//            p1, // Вес
//            bell // Функция влияния
//        );
//
//        fem_sol.save_as_vtk("heat.vtk", T);
//        save_raw_data(msh, T);
//
//        std::cout << fem_sol.integrate_solution(T) << std::endl;

    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return PetscFinalize();
}
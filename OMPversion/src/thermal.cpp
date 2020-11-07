#include "solvers/influence_functions.hpp"
#include "solvers/heat_equation_solver.hpp"

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(16);
        omp_set_num_threads(4);

        static constexpr double r = 0.1, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};
        nonlocal::heat::heat_equation_solver<double, int> fem_sol{msh};

        const auto T = fem_sol.stationary(
            { // Граничные условия
                {   // Down
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                },

                {   // Up
                    [](const std::array<double, 2>&) { return 1; },
                    nonlocal::heat::boundary_t::TEMPERATURE
                }
            }, 
            [](const std::array<double, 2>&) { return 0; }, // Правая часть
            r,  // Радиус влияния 
            p1, // Вес
            bell // Функция влияния
        );

        fem_sol.save_as_vtk("heat.vtk", T);

    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
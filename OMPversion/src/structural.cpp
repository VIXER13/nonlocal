#include "solvers/influence_functions.hpp"
#include "solvers/structural_solver.hpp"

int main(int argc, char** argv) {
    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(16);
        omp_set_num_threads(4);

        static constexpr double r = 0.05, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};

        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 42};
        nonlocal::structural::structural_solver<double, int> fem_sol{msh, params};

        const auto U = fem_sol.stationary(
            {
                {  // Down
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return -1; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {   // Up
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                },
            },
            r, p1, bell
        );

        const auto [strain, stress] = fem_sol.strains_and_stress(U, p1, bell);
        std::cout << "Energy U = " << fem_sol.calc_energy(strain, stress) << std::endl;

        fem_sol.save_as_vtk("structural.vtk", U, strain, stress);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
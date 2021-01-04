#include <iostream>
#include <petsc.h>
#include "influence_functions.hpp"
#include "structural_solver.hpp"

namespace {

template<class T>
void save_raw_data(const mesh::mesh_2d<T>& msh,
                   const std::vector<std::array<T, 3>>& strain,
                   const std::vector<std::array<T, 3>>& stress) {
    std::ofstream eps11{"eps11.csv"},
                  sigma11{"sigma11.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << strain[i][0] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << stress[i][0] << std::endl;
    }
}

}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(2);

        static constexpr double r = 0.2, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        mesh::mesh_2d<double> msh{argv[1]};

        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 21};
        nonlocal::structural::structural_solver<double, int> fem_sol{msh, params};

        const auto U = fem_sol.stationary(
            {
                {  // Down
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::DISPLACEMENT
                },

                {  // Right
                    [](const std::array<double, 2>&) { return 1; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {  // Up
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    nonlocal::structural::boundary_t::PRESSURE
                },

                {  // Left
                    [](const std::array<double, 2>&) { return 0; },
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    nonlocal::structural::boundary_t::PRESSURE
                },
            },
            r, p1, bell
        );

//        const auto U = fem_sol.stationary(
//            {
//                {  // Down
//                    [](const std::array<double, 2>&) { return 0; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    nonlocal::structural::boundary_t::PRESSURE
//                },
//
//                {  // Right
//                    [](const std::array<double, 2>&) { return 1; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    nonlocal::structural::boundary_t::PRESSURE
//                },
//
//                {  // Horizontal
//                    [](const std::array<double, 2>&) { return 0; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    nonlocal::structural::boundary_t::DISPLACEMENT
//                },
//
//                {  // Left
//                    [](const std::array<double, 2>&) { return -1; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    nonlocal::structural::boundary_t::PRESSURE
//                },
//
//                {  // Up
//                    [](const std::array<double, 2>&) { return 0; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    nonlocal::structural::boundary_t::PRESSURE
//                },
//
//                {  // Vertical
//                    [](const std::array<double, 2>&) { return 0; },
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::DISPLACEMENT,
//                    nonlocal::structural::boundary_t::PRESSURE
//                },
//            },
//            r, p1, bell
//        );

        const auto [strain, stress] = fem_sol.strains_and_stress(U, p1, bell);
        std::cout << "Energy U = " << fem_sol.calc_energy(strain, stress) << std::endl;

        save_raw_data(msh, strain, stress);
        fem_sol.save_as_vtk("structural.vtk", U, strain, stress);
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return PetscFinalize();
}
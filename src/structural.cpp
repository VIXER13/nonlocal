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
        std::cout.precision(7);

        static constexpr double r = 0.05, p1 = 1./2.;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_info = std::make_shared<mesh::mesh_info<double, int>>(mesh);
        if (p1 < 0.999)
            mesh_info->find_neighbours(1.2*r, mesh::balancing_t::SPEED);

        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 21};
        nonlocal::structural::structural_solver<double, int> fem_sol{mesh_info, params};


        // T
        const auto U = fem_sol.stationary(
            {
                {  // Down
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return -1; }
                },

                {  // Up
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; }
                }
            },
            p1, bell
        );


//        const auto U = fem_sol.stationary(
//            {
//                {  // Down
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return -1; }
//                },
//
//                {}, // Right
//
//                {   // Horizontal
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::DISPLACEMENT,
//                    [](const std::array<double, 2>&) { return 0; }
//                },
//
//                {}, // Left
//
//                {   // Vertical
//                    nonlocal::structural::boundary_t::DISPLACEMENT,
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 0; }
//                },
//
//                {   // Up
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 0; },
//                    nonlocal::structural::boundary_t::PRESSURE,
//                    [](const std::array<double, 2>&) { return 1; }
//                },
//
//                {} // Circle
//            },
//            p1, bell
//        );

        if (mesh_info->rank() == 0)
        {
            const auto [strain, stress] = fem_sol.strains_and_stress(U, p1, bell);
            std::cout << "Energy U = " << fem_sol.calc_energy(strain, stress) << std::endl;

            std::vector<double> strain11(strain.size()),
                                strain12(strain.size()),
                                strain22(strain.size());
            for(size_t i = 0; i < strain.size(); ++i) {
                strain11[i] = strain[i][0];
                strain12[i] = strain[i][1];
                strain22[i] = strain[i][2];
            }

            //save_raw_data(mesh, strain, stress);
            fem_sol.save_as_vtk("structural.vtk", U, strain, stress,
                                mesh_info->calc_gradient(strain11),
                                mesh_info->calc_gradient(strain12),
                                mesh_info->calc_gradient(strain22));
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    return PetscFinalize();
}
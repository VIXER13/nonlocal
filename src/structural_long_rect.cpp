#include <iostream>
//#include <petsc.h>
#include "influence_functions.hpp"
#include "structural_solver.hpp"

namespace {

template<class T>
void save_raw_data(const mesh::mesh_2d<T>& msh,
                   const nonlocal::structural::solution<T, int>& sol) {
    std::ofstream eps11{"eps11.csv"},
                  eps22{"eps22.csv"},
                  eps12{"eps12.csv"},
                  sigma11{"sigma11.csv"},
                  sigma22{"sigma22.csv"},
                  sigma12{"sigma12.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[0][i] << std::endl;
        eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[1][i] << std::endl;
        eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[2][i] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [0][i] << std::endl;
        sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [1][i] << std::endl;
        sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [2][i] << std::endl;
    }
}

double f1(const std::array<double, 2>& x) { return 1; }
double f2(const std::array<double, 2>& x) { return x[1] < 0.5 ? x[1] : 1 - x[1]; }
double f3(const std::array<double, 2>& x) { return f2(x) - 0.25; }
double f4(const std::array<double, 2>& x) { return 0.5 - f2(x); }
double f5(const std::array<double, 2>& x) { return 0.25 - f2(x); }
double f6(const std::array<double, 2>& x) { return x[1] < 0.5 ? 0.5 - x[1] : x[1] - 0.5; }
double f7(const std::array<double, 2>& x) { return x[1] < 0.45 || x[1] > 0.55 ? 0 : 10; }

}

int main(int argc, char** argv) {
    //PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);
    MPI_Init(&argc, &argv);

    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);

        static constexpr double r = 0.1, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<mesh::mesh_proxy<double, int>>(mesh);
        nonlocal::structural::structural_solver<double, int> fem_sol{mesh_proxy};
        if (p1 < 0.999)
            mesh_proxy->find_neighbours(1.3 * r, mesh::balancing_t::MEMORY);

        nonlocal::structural::calculation_parameters<double> params;
        params.nu = 0.3;
        params.E  = 21;

        auto sol = fem_sol.stationary(
            params,
            {
                {  // Right
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return f1(x); },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                {   // Horizontal
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Left
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return -f1(x); },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                {   // Vertical
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                }
            },
            {}, // Правая часть
            p1, bell
        );

        sol.calc_strain_and_stress();

        if (mesh_proxy->rank() == 0)
        {
            std::cout << "Energy = " << sol.calc_energy() << std::endl;
            sol.save_as_vtk("structural.vtk");
            save_raw_data(*mesh, sol);
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        return EXIT_FAILURE;
    }

    //return PetscFinalize();
    MPI_Finalize();
    return 0;
}
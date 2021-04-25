#include <iostream>
#include "influence_functions.hpp"
#include "structural_solver.hpp"

namespace {

template<class T>
void save_raw_data(const std::string& path,
                   const mesh::mesh_2d<T>& msh,
                   const nonlocal::structural::solution<T, int>& sol) {
    std::ofstream u1{path + "/u1.csv"},
            u2{path + "/u2.csv"},
            eps11{path + "/eps11.csv"},
            eps22{path + "/eps22.csv"},
            eps12{path + "/eps12.csv"},
            sigma11{path + "/sigma11.csv"},
            sigma22{path + "/sigma22.csv"},
            sigma12{path + "/sigma12.csv"};
    for(size_t i = 0; i < msh.nodes_count(); ++i) {
        u1      << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.displacement()[0][i] << std::endl;
        u2      << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.displacement()[1][i] << std::endl;
        eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[0][i] << std::endl;
        eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[1][i] << std::endl;
        eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.strains()[2][i] << std::endl;
        sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [0][i] << std::endl;
        sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [1][i] << std::endl;
        sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << sol.stress() [2][i] << std::endl;
    }
}

}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 5) {
        std::cerr << "Input format [program name] <path to mesh> <r> <p1> <output_path>";
        PetscFinalize();
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);

        const double r = std::stod(argv[2]), p1 = std::stod(argv[3]);
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);

        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_proxy = std::make_shared<mesh::mesh_proxy<double, int>>(mesh);
        if (p1 < 0.999)
            mesh_proxy->find_neighbours(r, mesh::balancing_t::MEMORY);

        nonlocal::structural::structural_solver<double, int, long long> fem_sol{mesh_proxy};
        nonlocal::structural::calculation_parameters<double> parameters;
        parameters.nu = 0.3;
        parameters.E = 21;

        auto sol = fem_sol.stationary(parameters,
              { // Граничные условия
                      {  // Up
                              nonlocal::structural::boundary_t::DISPLACEMENT,
                              [](const std::array<double, 2>&) { return 0; },
                              nonlocal::structural::boundary_t::DISPLACEMENT,
                              [](const std::array<double, 2>&) { return 0; }
                      },

                      {   // Down
                              nonlocal::structural::boundary_t::PRESSURE,
                              [](const std::array<double, 2>&) { return 0; },
                              nonlocal::structural::boundary_t::PRESSURE,
                              [](const std::array<double, 2>&) { return -1; }
                      }
              },
              {}, // Правая часть
              p1, // Вес
              bell // Функция влияния
        );

        sol.calc_strain_and_stress();

        if (mesh_proxy->rank() == 0) {
            std::cout << "Energy = " << sol.calc_energy() << std::endl;
            using namespace std::string_literals;
            sol.save_as_vtk(argv[4] + "/structural.vtk"s);
            save_raw_data(argv[4], *mesh, sol);
        }
    } catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
        PetscFinalize();
        return EXIT_FAILURE;
    } catch(...) {
        std::cerr << "Unknown error." << std::endl;
        PetscFinalize();
        return EXIT_FAILURE;
    }

    PetscFinalize();
    return EXIT_SUCCESS;
}
#include <iostream>
#include <petsc.h>
#include "influence_functions.hpp"
#include "heat_equation_solver.hpp"
#include "structural_solver.hpp"

namespace {

    template<class T>
    void save_raw_data(const mesh::mesh_2d<T>& msh,
                       const std::vector<std::array<T, 3>>& strain,
                       const std::vector<std::array<T, 3>>& stress) {
        std::ofstream eps11{"eps11.csv"},
                eps22{"eps22.csv"},
                eps12{"eps12.csv"},
                sigma11{"sigma11.csv"},
                sigma22{"sigma22.csv"},
                sigma12{"sigma12.csv"};
        for(size_t i = 0; i < msh.nodes_count(); ++i) {
            eps11   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << strain[i][0] << std::endl;
            eps22   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << strain[i][1] << std::endl;
            eps12   << msh.node(i)[0] << "," << msh.node(i)[1] << "," << strain[i][2] << std::endl;
            sigma11 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << stress[i][0] << std::endl;
            sigma22 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << stress[i][1] << std::endl;
            sigma12 << msh.node(i)[0] << "," << msh.node(i)[1] << "," << stress[i][2] << std::endl;
        }
    }

    double f1(const std::array<double, 2>& x) { return 1; }
    double f2(const std::array<double, 2>& x) { return x[1] < 0.5 ? x[1] : 1 - x[1]; }
    double f3(const std::array<double, 2>& x) { return f2(x) - 0.25; }
    double f4(const std::array<double, 2>& x) { return 0.5 - f2(x); }
    double f5(const std::array<double, 2>& x) { return 0.25 - f2(x); }
    double f6(const std::array<double, 2>& x) { return x[1] < 0.5 ? 0.5 - x[1] : x[1] - 0.5; }

}

int main(int argc, char** argv) {
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr); CHKERRQ(ierr);

    if(argc < 2) {
        std::cerr << "Input format [program name] <path to mesh>";
        return EXIT_FAILURE;
    }

    try {
        std::cout.precision(7);

        static constexpr double r = 0.1, p1 = 0.5;
        static const nonlocal::influence::polynomial<double, 2, 1> bell(r);
        auto mesh = std::make_shared<mesh::mesh_2d<double>>(argv[1]);
        auto mesh_info = std::make_shared<mesh::mesh_info<double, int>>(mesh);
        if (p1 < 0.999)
            mesh_info->find_neighbours(1.3 * r, mesh::balancing_t::SPEED);

        const nonlocal::structural::parameters<double> params = {.nu = 0.3, .E = 21};
        nonlocal::heat::heat_equation_solver<double, int> fem_sol_thermal{mesh_info};
        nonlocal::structural::structural_solver<double, int> fem_sol{mesh_info, params};

        auto T = fem_sol_thermal.stationary(
            { // Граничные условия
                { // Down
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                { // Left
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                { // Hor
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                { // Right
                    nonlocal::heat::boundary_t::TEMPERATURE,
                    [](const std::array<double, 2>& x) { return 1; },
                },

                { // Up
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },

                { // Ver
                    nonlocal::heat::boundary_t::FLOW,
                    [](const std::array<double, 2>& x) { return 0; },
                },
            },
            {[](const std::array<double, 2>& x) { return 0; }}, // Правая часть
            p1, bell
        );

        if(mesh_info->rank() == 0) {
            T.save_as_vtk("heat.vtk");
        }

        auto sol = fem_sol.stationary(
            {
                { // Down
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Left
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::DISPLACEMENT,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Hor
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Right
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Up
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>& x) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                },

                { // Ver
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; },
                    nonlocal::structural::boundary_t::PRESSURE,
                    [](const std::array<double, 2>&) { return 0; }
                }
            },
            {}, // Правая часть
            0.01, T.get_temperature(),
            p1, bell
        );

        if (mesh_info->rank() == 0)
        {
            sol.calc_strain_and_stress();
            sol.save_as_vtk("structural.vtk");
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
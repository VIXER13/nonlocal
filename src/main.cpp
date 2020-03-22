#include <iostream>
#include "nonlocal_influence_functions.hpp"
#include "mesh_2d.hpp"
#include "heat_equation_solver.hpp"
#include "static_equation_solver.hpp"
#include "thermomechanical_equation_solver.hpp"
#include "Eigen/Core"
#include "omp.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

using heat_boundary_type = heat_equation_with_nonloc::boundary_type;
using mechanical_boundary_type = static_equation_with_nonloc::boundary_type;

int main()
{
    std::cout.precision(5);
    const int threads = 4;
    omp_set_num_threads(threads);
    Eigen::initParallel();
    Eigen::setNbThreads(threads);

    const double r = 0.1, p1 = 1;
    //influence_function::constant<double> constant_fun(r);
    influence_function::polinomial<double, 1, 1> bell11(r);
    //influence_function::normal_distribution<double> norm(r);

    mesh_2d<double> mesh(mesh_2d<double>::BILINEAR, 100, 100, 1., 1.);
    mesh.find_neighbors_for_elements(1.05*r);

    size_t neighbors_count = 0;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        neighbors_count += mesh.neighbor(i).size();
    std::cout << "Average number of neighbors: " << double(neighbors_count) / mesh.elements_count() << std::endl;

    const auto T = heat_equation_with_nonloc::stationary<double, int>(
        mesh,
        { { // DOWN
            [](double x, double y) { return 0; },
            heat_boundary_type::FLOW },

          { // RIGHT
            [](double x, double y) { return 200; },
            heat_boundary_type::TEMPERATURE }, 

          { // UP
            [](double x, double y) { return 0; },
            heat_boundary_type::FLOW }, 

          { // LEFT
            [](double x, double y) { return 0; },
            heat_boundary_type::TEMPERATURE } },
        [](double, double) { return 0; },
        p1, bell11, 0.);
  
    const auto u = thermomechanical_equation_with_nonloc::stationary<double, int>(
        mesh, {.nu = 0.3, .E = 2.1e5, .alpha = 13e-6},
        { { // DOWN
            [](double, double) { return 0; },
            [](double, double) { return 0; },
            mechanical_boundary_type::PRESSURE,
            mechanical_boundary_type::PRESSURE },

          { // RIGHT
            [](double, double) { return 0; },
            [](double, double) { return 0; },
            mechanical_boundary_type::PRESSURE,
            mechanical_boundary_type::PRESSURE },

          { // UP
            [](double, double) { return 0; },
            [](double, double) { return 0; },
            mechanical_boundary_type::PRESSURE,
            mechanical_boundary_type::PRESSURE },

          { // LEFT
            [](double, double) { return 0; },
            [](double, double) { return 0; },
            mechanical_boundary_type::DISPLACEMENT,
            mechanical_boundary_type::DISPLACEMENT } },

        { [](double, double) { return 0; },
          [](double, double) { return 0; } },
          T, p1, bell11);

    mesh.find_neighbors_for_nodes(1.05*r);
    auto [eps11, eps22, eps12, sigma11, sigma22, sigma12] = 
        static_equation_with_nonloc::strains_and_stress<double, int>(mesh, u, {.nu = 0.3, .E = 2.1e5}, p1, bell11);
    static_equation_with_nonloc::raw_output("results//", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);
    static_equation_with_nonloc::save_as_vtk("results//loc.vtk", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);

    return 0;
}
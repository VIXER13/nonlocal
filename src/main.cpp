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

int main(int argc, char** argv)
{
    std::cout.precision(5);
    const int threads = 4;
    omp_set_num_threads(threads);
    Eigen::initParallel();
    Eigen::setNbThreads(threads);

    const double r = 10, p1 = 0.66;
    //nonlocal::influence::polynomial<double, 2, 1> bell11(r);
    nonlocal::influence::normal_distribution bell11(r);

    //mesh_2d<double> mesh(mesh_2d<double>::BILINEAR, 50, 50, 1., 1.);
    mesh_2d<double> mesh(argv[1]);
    mesh.find_neighbors_for_elements(3.*r);

    size_t neighbors_count = 0;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        neighbors_count += mesh.neighbor(i).size();
    std::cout << "Average number of neighbors: " << double(neighbors_count) / mesh.elements_count() << std::endl;
  
    const auto u = static_equation_with_nonloc::stationary<double, int>(
        mesh, {.nu = 0.2, .E = 1},
        { { // DOWN
                [](double, double) { return 0; },
                [](double, double) { return 0; },
                mechanical_boundary_type::DISPLACEMENT,
                mechanical_boundary_type::DISPLACEMENT },

            { // UP
                [](double, double) { return 0; },
                [](double, double) { return 1; },
                mechanical_boundary_type::PRESSURE,
                mechanical_boundary_type::DISPLACEMENT },
        },
        
        { [](double, double) { return 0; },
          [](double, double) { return 0; } },

        p1, bell11);

    mesh.find_neighbors_for_nodes(3.*r);
    auto [eps11, eps22, eps12, sigma11, sigma22, sigma12] = 
        static_equation_with_nonloc::strains_and_stress<double, int>(mesh, u, {.nu = 0.2, .E = 1}, p1, bell11);
    static_equation_with_nonloc::raw_output("results//", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);
    static_equation_with_nonloc::save_as_vtk("results//loc.vtk", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);

    return 0;
}
#include <iostream>
#include "nonlocal_influence_functions.hpp"
#include "mesh_2d.hpp"
#include "heat_equation_solver.hpp"
#include "static_analysis.hpp"
#include "Eigen/Core"
#include "omp.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

int main()
{
    std::cout.precision(5);
    const int threads = 4;
    omp_set_num_threads(threads);
    Eigen::initParallel();
    Eigen::setNbThreads(threads);

    const double r = 0.1;
    influence_function::constant<double> constant_fun(r);
    influence_function::polinomial<double, 1, 1> bell11(r);
    influence_function::normal_distribution<double> norm(r);

    mesh_2d<double> mesh(mesh_2d<double>::BILINEAR, 5, 5, 1., 1.);
    mesh.find_neighbors_for_elements(1.05*r);

    size_t neighbors_count = 0;
    for(size_t i = 0; i < mesh.elements_count(); ++i)
        neighbors_count += mesh.neighbor(i).size();
    std::cout << "Average number of neighbors: " << double(neighbors_count) / mesh.elements_count() << std::endl;
    
    
    {
    using namespace statics_with_nonloc;
    const parameters<double> param = {.nu = 0.3, .E = 2.1e5};
    const double p1 = 1.;
    
    Eigen::VectorXd u = stationary(mesh, param,
                                   { { [](double, double) { return 0; },
                                       [](double, double) { return 0.; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },

                                     { [](double, double) { return 5000.; },
                                       [](double, double) { return 0.; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },

                                     { [](double, double) { return 0; },
                                       [](double, double) { return 0; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },

                                     { [](double, double) { return 0.; },
                                       [](double, double) { return 0.; },
                                       boundary_type::DISPLACEMENT,
                                       boundary_type::DISPLACEMENT } },
                                   p1, bell11);
    
    /*
    Eigen::VectorXd u = stationary(mesh, param,
                                   { { [](double, double) { return 0.; },
                                       [](double, double) { return 0.; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },

                                     { [](double, double) { return 5000; },
                                       [](double, double) { return 0.; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },

                                     { [](double, double) { return 0; },
                                       [](double, double) { return 0; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE },
                                    
                                     { [](double, double) { return -5000.; },
                                       [](double, double) { return 0.; },
                                       boundary_type::PRESSURE,
                                       boundary_type::PRESSURE, } },
                                   p1, bell11);
                                   */

    mesh.find_neighbors_for_nodes(1.05*r);
    auto [eps11, eps22, eps12, sigma11, sigma22, sigma12] = strains_and_stress(mesh, u, param, p1, bell11);
    raw_output("results//", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);
    save_as_vtk("results//loc.vtk", mesh, u, eps11, eps22, eps12, sigma11, sigma22, sigma12);
    }
    

    std::cout << std::endl << std::endl;
    /*
    {
    using namespace heat_equation_with_nonloc;
    stationary<double, int>(std::string("results//test.csv"), mesh,
                        { { [](double x, double y) { return x*x+y*y; }, boundary_type::TEMPERATURE }, 
                          { [](double x, double y) { return x*x+y*y; }, boundary_type::TEMPERATURE      }, 
                          { [](double x, double y) { return x*x+y*y; }, boundary_type::TEMPERATURE }, 
                          { [](double x, double y) { return x*x+y*y; }, boundary_type::TEMPERATURE      } },
                        [](double, double) { return -4.; },
                        1., bell11, 0.);
    }
    */

    /*
    heat_equation_with_nonloc::nonstationary(std::string("results//nonstationary_test//"), mesh, 0.01, 100,
                                            { { boundary_type::TEMPERATURE, [](double, double) { return  0.; } }, 
                                              { boundary_type::FLOW,        [](double, double) { return  0.; } }, 
                                              { boundary_type::TEMPERATURE, [](double, double) { return  1.; } }, 
                                              { boundary_type::FLOW,        [](double, double) { return  0.; } } },
                                            [](double, double) { return 0.; },
                                            [](double, double) { return 0.; },
                                            .5, bell11, 1);
    */

    return 0;
}
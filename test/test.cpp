#include <iostream>
#include "solvers/heat_equation_solver.hpp"

int main(int argc, char** argv) {
    mesh::mesh_2d<double> mesh{argv[1]};

    std::cout << "nodes count: " << mesh.nodes_count() << std::endl;
    for(size_t i = 0; i < mesh.nodes_count(); ++i)
        std::cout << i << " " << mesh.node(i)[0] << " " << mesh.node(i)[1] << std::endl;
    std::cout << std::endl;

    std::cout << "elements count: " << mesh.elements_count() << std::endl;
    for(size_t el = 0; el < mesh.elements_count(); ++el) {
        const auto& e = mesh.element_2d(mesh.element_2d_type(el));
        std::cout << el;
        for(size_t i = 0; i < e->nodes_count(); ++i)
            std::cout << " " << mesh.node_number(el, i);
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "boundary groups count: " << mesh.boundary_groups_count() << std::endl;
    for(size_t b = 0; b < mesh.boundary_groups_count(); ++b) {
        std::cout << b << " bound " << mesh.elements_count(b) << std::endl;
        for(size_t el = 0; el < mesh.elements_count(b); ++el) {
            const auto& e = mesh.element_1d(mesh.element_1d_type(b, el));
            for(size_t i = 0; i < e->nodes_count(); ++i)
                std::cout << " " << mesh.node_number(b, el, i);
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    // const auto shifts = nonlocal::quadrature_shifts_init(mesh);
    // for(size_t i = 0; i < shifts.size(); ++i)
    //     std::cout << shifts[i] << ' ';
    
    // const auto coords = nonlocal::approx_all_quad_nodes(mesh, shifts);
    // for(size_t i = 0; i < coords.size(); ++i)
    //     std::cout << coords[i][0] << ' ' << coords[i][1] << std::endl;

    // std::cout << std::endl;

    // const auto jacobi = nonlocal::approx_all_jacobi_matrices(mesh, shifts);
    // for(size_t i = 0; i < jacobi.size(); ++i)
    //     std::cout << jacobi[i][0] << ' ' << jacobi[i][1] << ' ' << jacobi[i][2] << ' ' << jacobi[i][3] << std::endl;

    return 0;
}
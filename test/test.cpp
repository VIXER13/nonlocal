#include <iostream>
#include "include/mesh/mesh.hpp"
#include "include/containers/matrix.hpp"

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
        std::cout << b << " bound " << mesh.nodes_count(b) << std::endl;
        for(size_t el = 0; el < mesh.nodes_count(b); ++el) {
            const auto& e = mesh.element_1d(mesh.element_1d_type(b, el));
            for(size_t i = 0; i < e->nodes_count(); ++i)
                std::cout << " " << mesh.node_number(b, el, i);
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    const auto& e = mesh.element_2d(mesh.element_2d_type(0));
    matrix<double> jacobi_matrices(e->qnodes_count(), 4, 0);
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            jacobi_matrices(q, 0) += mesh.node(mesh.node_number(0, i))[0] * e->qNxi (i, q);
            jacobi_matrices(q, 1) += mesh.node(mesh.node_number(0, i))[0] * e->qNeta(i, q);
            jacobi_matrices(q, 2) += mesh.node(mesh.node_number(0, i))[1] * e->qNxi (i, q);
            jacobi_matrices(q, 3) += mesh.node(mesh.node_number(0, i))[1] * e->qNeta(i, q);
        }
    std::cout << jacobi_matrices << std::endl;

    return 0;
}
#include <iostream>
#include "su2_parser.hpp"

namespace {

using T = double;
using I = int64_t;

}

int main(int argc, char** argv) {
    nonlocal::mesh::mesh_container_2d<double, int64_t> mesh{argv[1]};

    std::cout << "groups_names_1d: ";
    for(const std::string& name : mesh.groups_names_1d())
        std::cout << name << ' ';
    std::cout << std::endl;

    std::cout << "groups_names_2d: ";
    for(const std::string& name : mesh.groups_names_2d())
        std::cout << name << ' ';
    std::cout << std::endl;

    std::cout << "elements_count = " << mesh.elements_count() << std::endl;
    std::cout << "elements_1d_count = " << mesh.elements_1d_count() << std::endl;
    std::cout << "elements_2d_count = " << mesh.elements_2d_count() << std::endl;

    std::cout << "nodes_count = " << mesh.nodes_count() << std::endl;
    // for(const size_t node : std::ranges::iota_view{0u, mesh.nodes_count()}) {
    //     const auto& coord = mesh.node_coord(node);
    //     std::cout << "node = " << node << " x = " << coord[0] << " y = " << coord[1] << std::endl;
    // }

    return 0;
}
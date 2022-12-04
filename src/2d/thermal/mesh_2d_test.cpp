#include <iostream>
#include "mesh_2d.hpp"
#include "thermal_conductivity_matrix_2d.hpp"
#include "influence_functions_2d.hpp"

namespace {

using T = double;
using I = int64_t;

void test_mesh(const nonlocal::mesh::mesh_container_2d<T, I>& mesh) {
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
}

}

int main(int argc, char** argv) {
    std::cout.precision(2);
    auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
    test_mesh(mesh->container());

    //const auto process_nodes = mesh->process_nodes();
    //std::cout << process_nodes.front() << " " << process_nodes.back() << std::endl;

    //nonlocal::thermal::equation_parameters_2d<T, nonlocal::material_t::ISOTROPIC> param;
    nonlocal::thermal::thermal_conductivity_matrix_2d<T, I, I> matrix{mesh};
    matrix.compute({1.}, nonlocal::material_t::ISOTROPIC, std::vector<bool>(mesh->container().nodes_count(), true), 1., nonlocal::influence::constant_2d<T>{1}, true);
    std::cout << Eigen::MatrixXd{matrix.matrix_inner()} << std::endl;


    return 0;
}
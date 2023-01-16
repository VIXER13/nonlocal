#include <iostream>
//#include "stationary_heat_equation_solver_2d.hpp"
#include <mesh_2d.hpp>
#include <thermal_boundary_conditions_2d.hpp>
#include "base/solvers_utils.hpp"
#include "influence_functions_2d.hpp"

namespace {

using T = double;
using I = int64_t;

void test_mesh(const nonlocal::mesh::mesh_container_2d<T, I>& mesh) {
    std::cout << "groups_names_1d: ";
    for(const std::string& name : mesh.groups_1d())
        std::cout << name << ' ';
    std::cout << std::endl;

    std::cout << "groups_names_2d: ";
    for(const std::string& name : mesh.groups_2d())
        std::cout << name << ' ';
    std::cout << std::endl;
    
    std::cout << "elements_1d_count = " << mesh.elements_1d_count() << std::endl;
    std::cout << "elements_2d_count = " << mesh.elements_2d_count() << std::endl;

    std::cout << "nodes_count = " << mesh.nodes_count() << std::endl;
    // for(const size_t node : std::ranges::iota_view{0u, mesh.nodes_count()}) {
    //     const auto& coord = mesh.node_coord(node);
    //     std::cout << "node = " << node << " x = " << coord[0] << " y = " << coord[1] << std::endl;
    // }
}

template<class U, nonlocal::physics_t Physics, size_t DoF>
void foo(const nonlocal::boundaries_conditions_2d<U, Physics, DoF>&) {
//void foo(const std::unordered_map<std::string, std::unique_ptr<nonlocal::boundary_condition_2d<U, Physics>>>&) {
    std::cout << "DoF = " << DoF << std::endl;
    std::cout << "Phys = " << int(Physics) << std::endl;
    //std::cout << int(Physics) << std::endl;
}

}

int main(const int argc, const char *const *const argv) {
    std::cout.precision(2);
    auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
    test_mesh(mesh->container());

    nonlocal::boundaries_conditions_2d<T, nonlocal::physics_t::THERMAL, 1> boundary_conditions;
    //std::unordered_map<std::string, std::unique_ptr<nonlocal::boundary_condition_2d<T, nonlocal::physics_t::THERMAL>>> boundary_conditions;
    boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(-1);
    boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(-1);

    // foo(boundary_conditions);

    std::cout << "before inner_nodes" << std::endl;
    auto inner_nodes = nonlocal::utils::inner_nodes(mesh->container(), boundary_conditions);
    std::cout << "after inner_nodes" << std::endl;

    for(size_t i : std::ranges::iota_view{0u, inner_nodes.size()})
        std::cout << i << ' ' << (inner_nodes[i] ? "true" : "false") << std::endl;

    // std::cout << "before boundaries" << std::endl;

    // nonlocal::thermal::boundaries_conditions_2d<T> boundary_conditions;
    // boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::convection_2d<T>>(
    //     -1., 1.
    // );
    // boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
    //     1.
    // );
    // // boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    // //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // // );
    // // boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    // //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // // );
    // // boundary_conditions["Up"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    // //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // // );
    // // boundary_conditions["Down"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    // //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // // );

    // static constexpr auto right_part = [](const std::array<T, 2>& x) {
    //     return 0;
    // };

    // std::cout << "before tast" << std::endl;
    // nonlocal::thermal::parameter_2d<T> parameters;
    // auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
    //     mesh, parameters, boundary_conditions, right_part, 1., nonlocal::influence::constant_2d<T>{1}
    // );

    // solution.calc_flux();
    // solution.save_as_vtk("./test.vtk");
    // nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./T.csv"}, mesh->container(), solution.temperature());
    // nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./TX.csv"}, mesh->container(), solution.flux()[0]);
    // nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./TY.csv"}, mesh->container(), solution.flux()[1]);

    return 0;
}
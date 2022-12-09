#include <iostream>
#include "stationary_heat_equation_solver_2d.hpp"
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

}

int main(const int argc, const char *const *const argv) {
    std::cout.precision(2);
    auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<T, I>>(argv[1]);
    test_mesh(mesh->container());


    std::cout << "before boundaries" << std::endl;

    nonlocal::thermal::boundaries_conditions_2d<T> boundary_conditions;
    boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        -1.
    );
    boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::flux_2d<T>>(
        1.
    );
    // boundary_conditions["Left"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // );
    // boundary_conditions["Right"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // );
    // boundary_conditions["Up"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // );
    // boundary_conditions["Down"] = std::make_unique<nonlocal::thermal::temperature_2d<T>>(
    //     [](const std::array<T, 2>& x) constexpr noexcept { return x[0] * x[0] + x[1] * x[1]; }
    // );

    static constexpr auto right_part = [](const std::array<T, 2>& x) {
        return -4;
    };

    std::cout << "before tast" << std::endl;
    nonlocal::thermal::parameter_2d<T> parameters;
    auto solution = nonlocal::thermal::stationary_heat_equation_solver_2d<I>(
        mesh, parameters, boundary_conditions, right_part, 1., nonlocal::influence::constant_2d<T>{1}
    );

    solution.calc_flux();
    solution.save_as_vtk("./test.vtk");
    nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./T.csv"}, mesh->container(), solution.temperature());
    nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./TX.csv"}, mesh->container(), solution.flux()[0]);
    nonlocal::mesh::utils::save_as_csv(std::filesystem::path{"./TY.csv"}, mesh->container(), solution.flux()[1]);

    return 0;
}
#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <mesh/mesh_2d/search_function.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

#include <embedded_files/solid_ring_su2.h>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal;
using namespace unit_tests;
using namespace mesh;
using namespace solver_2d;
using namespace nonlocal::solver_2d::influence;
using namespace mechanical;

constexpr T Expected_Error = T{0};
constexpr T Radius = T{0.2};

// TODO: Add analytical solution for the isotropic Lame solid ring problem
    
const suite<"nonlocal_example"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    mesh->neighbours(
        mesh::find_neighbours(*mesh, 
            {{"DEFAULT", {powered_distance<T>{Radius}, {Radius, Radius}}}}
        )
    );
    const raw_mechanical_parameters<T> parameters = { 
        {"DEFAULT", { 
            .model = { .influence = fast_polynomial<T, powered_distance<T>>{{Radius, Radius}}, .local_weight = 1 },
            .physical = { .elastic = isotropic_elastic_parameters<T>{ .young_modulus = 350., .poissons_ratio = 0.25 } } 
        }}
    };
    mechanical_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Horizontal"] = {
        nullptr,
        std::make_unique<displacement_2d<T>>(T{0})
    };
    boundaries_conditions["Vertical"] = {
        std::make_unique<displacement_2d<T>>(T{0}),
        nullptr
    };
    static constexpr T Inner_Pressure = T{0.05};
    boundaries_conditions["Inner"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return Inner_Pressure * point[X] / std::hypot(point[X], point[Y]); }),
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return Inner_Pressure * point[Y] / std::hypot(point[X], point[Y]); })
    };
    static constexpr T Outer_Pressure = T{0.01};
    boundaries_conditions["Outer"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return -Outer_Pressure * point[X] / std::hypot(point[X], point[Y]); }),
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return -Outer_Pressure * point[Y] / std::hypot(point[X], point[Y]); })
    };
    const auto solution = equilibrium_equation(mesh, parameters, boundaries_conditions);
};

}
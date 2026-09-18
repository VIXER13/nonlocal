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
constexpr T Inner_Pressure = T{0.05};
constexpr T Outer_Pressure = T{0.01};
constexpr T Inner_Radius = T{0.5};
constexpr T Outer_Radius = T{1.0};
constexpr T Er = T{50};
constexpr T Ef = T{200};
constexpr T nu_rf = T{0.05};
constexpr T nu_fr = Ef * nu_rf / Er;

constexpr T Radius = T{0.2};
    
const suite<"nonlocal_test"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    mesh->neighbours(
        mesh::find_neighbours(*mesh, 
            {{"DEFAULT", {powered_distance<T>{Radius}, {Radius, Radius}}}}
        )
    );

    const anisotropic_elastic_parameters<T> elastic = {
        .main_parameters = {
            .young_modulus = {Er, Ef},
            .poissons_ratio = {nu_rf, nu_fr},
            .shear_modulus = 180.
        },
        .angle = [](const std::array<T, 2>& point) noexcept { return std::atan2(point[Y], point[X]); }
    };
    const raw_mechanical_parameters<T> parameters = { { "DEFAULT", {
        .model = { .influence = fast_polynomial<T, powered_distance<T>>{{Radius, Radius}}, .local_weight = 0.5 },
        .physical = { .elastic = elastic }
    } } };
    mechanical_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Horizontal"] = {
        nullptr,
        std::make_unique<displacement_2d<T>>(T{0})
    };
    boundaries_conditions["Vertical"] = {
        std::make_unique<displacement_2d<T>>(T{0}),
        nullptr
    };
    boundaries_conditions["Inner"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return Inner_Pressure * std::cos(std::atan2(point[Y], point[X])); }),
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return Inner_Pressure * std::sin(std::atan2(point[Y], point[X])); })
    };
    boundaries_conditions["Outer"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return -Outer_Pressure * std::cos(std::atan2(point[Y], point[X])); }),
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return -Outer_Pressure * std::sin(std::atan2(point[Y], point[X])); })
    };
    std::cerr << "with preconditioner" << std::endl;
    const auto solution = equilibrium_equation(mesh, parameters, boundaries_conditions, std::vector<T>{}, std::function<std::array<T, 2>(const std::array<T, 2>&)>{}, true);
    std::cerr << "without preconditioner" << std::endl;
    const auto solution2 = equilibrium_equation(mesh, parameters, boundaries_conditions, std::vector<T>{}, std::function<std::array<T, 2>(const std::array<T, 2>&)>{}, false);
};

}
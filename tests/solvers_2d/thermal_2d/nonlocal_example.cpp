#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <mesh/mesh_2d/search_function.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>
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
using namespace solver_2d::thermal;

constexpr T Expected_Error = T{0};
constexpr T Inner_Temperature = T{0};
constexpr T Outer_Temperature = T{1};
constexpr T Norm_Temperature = Outer_Temperature / (Outer_Temperature - Inner_Temperature);
constexpr T Inner_Radius = T{0.5};
constexpr T Outer_Radius = T{1};
constexpr T Norm_Radius = Inner_Radius / Outer_Radius;
const T Coeff = T{1} / std::log(T{1} / Norm_Radius);
constexpr T Radius = T{0.2};

const suite<"thermal_isotropic_solid_ring"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    mesh->neighbours(
        mesh::find_neighbours(*mesh, 
            {{"DEFAULT", {powered_distance<T>{Radius}, {Radius, Radius}}}}
        )
    );
    const raw_thermal_parameters<T> parameters = {{"DEFAULT", {
        .model = { .influence = fast_polynomial<T, powered_distance<T>>{{Radius, Radius}}, .local_weight = 0.5 },
        .physical = {.conductivity = [](const std::array<T, 2>& point) { return T{1}; }}
    }}};
    thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Inner"] = std::make_unique<temperature_2d<T>>([](const std::array<T, 2>& point) { return Inner_Temperature; });
    boundaries_conditions["Outer"] = std::make_unique<temperature_2d<T>>([](const std::array<T, 2>& point) { return Outer_Temperature; });
    const auto solution = stationary_heat_equation_solver_2d(mesh, parameters, boundaries_conditions, {});
};

}
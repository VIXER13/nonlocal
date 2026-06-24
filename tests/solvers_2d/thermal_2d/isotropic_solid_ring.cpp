#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

#include <embedded_files/solid_ring_su2.h>

namespace {

using T = double;
using I = int64_t;
using namespace boost::ut;
using namespace nonlocal;
using namespace unit_tests;
using namespace mesh;
using namespace solver_2d::thermal;

constexpr T Expected_Error = T{0};
constexpr T Inner_Temperature = T{0};
constexpr T Outer_Temperature = T{1};
constexpr T Norm_Temperature = Outer_Temperature / (Outer_Temperature - Inner_Temperature);
constexpr T Inner_Radius = T{0.5};
constexpr T Outer_Radius = T{1};
constexpr T Norm_Radius = Inner_Radius / Outer_Radius;
const T Coeff = T{1} / std::log(T{1} / Norm_Radius);

const suite<"thermal_isotropic_solid_ring"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T, I>>(stream, mesh_format::SU2);
    const parameters_2d<T> parameters = {{"DEFAULT", {.physical = {.conductivity = T{1}}}}};
    thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Inner"] = std::make_unique<temperature_2d<T>>([](const std::array<T, 2>& point) { return Inner_Temperature; });
    boundaries_conditions["Outer"] = std::make_unique<temperature_2d<T>>([](const std::array<T, 2>& point) { return Outer_Temperature; });
    const auto solution = stationary_heat_equation_solver_2d<I>(mesh, parameters, boundaries_conditions, {});

    "temperature"_test = [&mesh, &solution] {
        static constexpr auto Expected_Temperature = [](const std::array<T, 2>& point) {
            return Norm_Temperature + Coeff * std::log(std::hypot(point[X], point[Y]));
        };
        static constexpr T Epsilon = 1.4e-3;
        const T error = norm_error(solution.temperature(), mesh->container(), Expected_Temperature);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "flux_x"_test = [&mesh, &solution] {
        static constexpr auto Expected_Flux_X = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            return -x * Coeff / (x * x + y * y);
        };
        static constexpr T Epsilon = 1.7e-2;
        const T error = norm_error(solution.flux()[X], mesh->container(), Expected_Flux_X);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "flux_y"_test = [&mesh, &solution] {
        static constexpr auto Expected_Flux_Y = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            return -y * Coeff / (x * x + y * y);
        };
        static constexpr T Epsilon = 1.7e-2;
        const T error = norm_error(solution.flux()[Y], mesh->container(), Expected_Flux_Y);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
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
using namespace solver_2d::thermal;
using namespace metamath::constants;
using namespace metamath::functions;

constexpr T Expected_Error = T{0};
constexpr T Inner_Radius = T{0.5};
constexpr T Outer_Radius = T{1};
constexpr T A = T{100};
constexpr T B = T{-50};
constexpr T Inner_Temperature = T{100};
constexpr T Emissivity = T{0.8};
constexpr T Heat_Transfer = T{0};
constexpr T Ambient_Temperature = T{0};

constexpr T temperature(const std::array<T, 2>& x) noexcept {
    const T r = std::hypot(x[X], x[Y]) - Inner_Radius;
    return Inner_Temperature + A * r + B * r * r;
}

constexpr T flux(const std::array<T, 2>& x) noexcept {
    const T r = std::hypot(x[X], x[Y]) - Inner_Radius;
    return -A - 2 * B * r;
}

const suite<"thermal_isotropic_solid_ring_radiation"> _ = [] {
    const T Outer_Temperature = temperature({T{0}, Outer_Radius});
    const T Outer_Flux = Emissivity * Stefan_Boltzmann_Constant<T> * power<4>(Outer_Temperature) - flux({T{0}, Outer_Radius});

    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    const parameters_2d<T> parameters = {{"DEFAULT", {.physical = {.conductivity = T{1}}}}};
    thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Inner"] = std::make_unique<temperature_2d<T>>(Inner_Temperature);
    boundaries_conditions["Outer"] = std::make_unique<combined_flux_2d<T>>(Outer_Flux, Heat_Transfer, Ambient_Temperature, Emissivity);
    const stationary_equation_parameters_2d<T> auxiliary_data = {
        .right_part = [](const std::array<T, 2>& point) { return -4 * B + (2 * B * Inner_Radius - A) / std::hypot(point[X], point[Y]); },
        .initial_distribution = [Outer_Temperature](const std::array<T, 2>& point) { return 0.5 * (Inner_Temperature + Outer_Temperature); },
        .max_iterations = 10
    };
    const auto solution = stationary_heat_equation_solver_2d(mesh, parameters, boundaries_conditions, auxiliary_data);

    "temperature"_test = [&mesh, &solution] {
        static constexpr T Epsilon = 4e-4;
        static constexpr auto Expected_Temperature = [](const std::array<T, 2>& point) { return temperature(point); };
        const T error = norm_error(solution.temperature(), mesh->container(), Expected_Temperature);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "flux"_test = [&mesh, &solution] {
        static constexpr auto Expected_Flux = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T coeff = flux(point) / std::hypot(x, y);
            return std::array{coeff * x, coeff * y};
        };
        static constexpr T Epsilon = 1.8e-2;
        const T error = norm_error(solution.flux(), mesh->container(), Expected_Flux);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
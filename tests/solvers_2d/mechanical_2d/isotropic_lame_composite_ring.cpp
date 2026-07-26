#include "coordinate_converter.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

#include <embedded_files/composite_ring_su2.h>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal;
using namespace unit_tests;
using namespace mesh;
using namespace solver_2d::mechanical;

constexpr T Expected_Error = T{0};
constexpr T Contact_Radius = T{0.75};
constexpr T Contact_Radius_Sqr = metamath::functions::power<2>(Contact_Radius);

// TODO: Add analytical solution for the composite ring problem

bool is_contact(const T r) noexcept {
    static constexpr T Epsilon = T{1e-12};
    return std::abs(r - Contact_Radius_Sqr) < Epsilon;
}

const suite<"isotropic_lame_composite_ring"> _ = [] {
    std::stringstream stream{ composite_ring_su2_data };
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    const raw_mechanical_parameters<T> parameters = {
        {"Inner_Material", {.physical = {.elastic = isotropic_elastic_parameters<T>{.young_modulus = 350., .poissons_ratio = 0.25 } } }},
        {"Outer_Material", {.physical = {.elastic = isotropic_elastic_parameters<T>{.young_modulus = 150., .poissons_ratio = 0.2  } } }}
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

    "displacement"_test = [&mesh, &solution] {
        static constexpr auto Expected_Displacement = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T r2 = x * x + y * y;
            const T inner = (T{2.5284} + r2) / (T{47303}  * r2);
            const T outer = (T{9.7572} - r2) / (T{140715} * r2);
            return r2 < Contact_Radius_Sqr ? std::array{x * inner, y * inner} : std::array{x * outer, y * outer};
        };
        static constexpr T Epsilon = 1.6e-3;
        const T error = norm_error(solution.displacement(), mesh->container(), Expected_Displacement);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "strain"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T r2 = x2 + y2;
            const T r4 = r2 * r2;
            const T mul = -x * y / r4;
            const T diff = (y2 - x2) / r4;
            using metamath::functions::power;
            using namespace metamath::operators;
            const auto Inner_Solution = [diff, mul]() { 
                return std::array{(T{1} + T{2.5284} * diff) / T{47303},
                                  (T{1} - T{2.5284} * diff) / T{47303},
                                                        mul / T{9354.31}};
            };
            const auto Outer_Solution = [diff, mul]() {
                return std::array{(T{-1} + T{9.7572} * diff) / T{140715},
                                  (T{-1} - T{9.7572} * diff) / T{140715},
                                                         mul / T{7210.83}};
            };
            if (is_contact(r2))
                return 0.5 * (Inner_Solution() + Outer_Solution());
            return r2 < Contact_Radius_Sqr ? Inner_Solution() : Outer_Solution();
        };
        static constexpr T Epsilon = 2.7e-2;
        const T error = norm_error(solution.strain(), mesh->container(), Expected_Strain);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T r2 = x2 + y2;
            const T r4 = r2 * r2;
            const T mul = -x * y / r4;
            const T diff = (y2 - x2) / r4;
            using metamath::functions::power;
            using namespace metamath::operators;
            const auto Inner_Solution = [diff, mul]() {
                return std::array{(T{1} + T{1.5170} * diff) / T{101.364},
                                  (T{1} - T{1.5170} * diff) / T{101.364},
                                                        mul / T{33.4082}};
            };
            const auto Outer_Solution = [diff, mul]() {
                return std::array{-(T{1} - T{6.5048} * diff) / T{750.481},
                                  -(T{1} + T{6.5048} * diff) / T{750.481},
                                                         mul / T{57.6866}};
            };
            if (is_contact(r2))
                return 0.5 * (Inner_Solution() + Outer_Solution());
            return r2 < Contact_Radius_Sqr ? Inner_Solution() : Outer_Solution();
        };
        static constexpr T Epsilon = 2.9e-2;
        const T error = norm_error(solution.stress(), mesh->container(), Expected_Stress);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
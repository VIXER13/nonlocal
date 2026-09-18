#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
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
using namespace solver_2d::mechanical;

constexpr T Expected_Error = T{0};

// TODO: Add analytical solution for the isotropic Lame solid ring problem
    
const suite<"isotropic_lame_solid_ring"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    const raw_mechanical_parameters<T> parameters = { 
        {"DEFAULT", { .physical = { .elastic = isotropic_elastic_parameters<T>{ .young_modulus = 350., .poissons_ratio = 0.25 } } }}
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
            const T hypot2 = x * x + y * y;
            const T value = (T{20} + T{3} * hypot2) / (T{420000} * hypot2);
            return std::array{x * value, y * value};
        };
        static constexpr T Epsilon = 1.1e-3;
        const T error = norm_error(solution.displacement(), mesh->container(), Expected_Displacement);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "strain"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            using metamath::functions::power;
            return std::array{
                (T{3} * power<2>(x2) + y2 * (T{ 20} + T{3} * y2) + x2 * (T{-20} + T{6} * y2)) / (T{420000} * power<2>(x2 + y2)),
                (T{3} * power<2>(x2) + y2 * (T{-20} + T{3} * y2) + x2 * (T{ 20} + T{6} * y2)) / (T{420000} * power<2>(x2 + y2)),
                -std::sin(2 * std::atan2(y, x)) / (T{21000} * (x * x + y * y))
            };
        };
        static constexpr T Epsilon = 3.1e-2;
        const T error = norm_error(solution.strain(), mesh->container(), Expected_Strain);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            using metamath::functions::power;
            return std::array{
                (power<2>(x2) + T{2} * x2 * (T{-2} + y2) + y2 * (T{ 4} + y2)) / (T{300} * power<2>(x2 + y2)),
                (power<2>(x2) + T{2} * x2 * (T{ 2} + y2) + y2 * (T{-4} + y2)) / (T{300} * power<2>(x2 + y2)),
                -std::sin(2 * std::atan2(y, x)) / (T{75} * (x * x + y * y))
            };
        };
        static constexpr T Epsilon = 3.3e-2;
        const T error = norm_error(solution.stress(), mesh->container(), Expected_Stress);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
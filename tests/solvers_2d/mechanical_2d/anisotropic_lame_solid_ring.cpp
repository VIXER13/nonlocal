#include "coordinate_converter.hpp"

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
constexpr T Inner_Pressure = T{0.05};
constexpr T Outer_Pressure = T{0.01};
constexpr T Inner_Radius = T{0.5};
constexpr T Outer_Radius = T{1.0};
constexpr T Er = T{50};
constexpr T Ef = T{200};
constexpr T nu_rf = T{0.05};
constexpr T nu_fr = Ef * nu_rf / Er;
    
const suite<"anisotropic_lame_solid_ring"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    const anisotropic_elastic_parameters<T> elastic = {
        .main_parameters = {
            .young_modulus = {Er, Ef},
            .poissons_ratio = {nu_rf, nu_fr},
            .shear_modulus = 180.
        },
        .angle = [](const std::array<T, 2>& point) noexcept { return std::atan2(point[Y], point[X]); }
    };
    const raw_mechanical_parameters<T> parameters = { { "DEFAULT", { .physical = { .elastic = elastic } } } };
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
    const auto solution = equilibrium_equation(mesh, parameters, boundaries_conditions);

    const T k = std::sqrt(Ef / Er);
    const T Delta = (T{1} - nu_rf * nu_fr) / (std::pow(Outer_Radius, 2 * k) - std::pow(Inner_Radius, 2 * k));
    const T A = Delta * (std::pow(Inner_Radius, k + 1) * Inner_Pressure - std::pow(Outer_Radius, k + 1) * Outer_Pressure) / (Er * (k + nu_fr));
    const T B = Delta * std::pow(Inner_Radius * Outer_Radius, k) * 
                (Inner_Radius * std::pow(Outer_Radius, k) * Inner_Pressure - Outer_Radius * std::pow(Inner_Radius, k) * Outer_Pressure) / (Er * (k - nu_fr));
    const auto Expected_Displacement_R = [A, B, k](const std::array<T, 2>& point) {
        const T r = std::hypot(point[X], point[Y]);
        return A * std::pow(r, k) + B * std::pow(r, -k);
    };

    const auto Expected_Polar_Strain = [A, B, k](const std::array<T, 2>& point) -> std::array<T, 3> {
        const T r = std::hypot(point[X], point[Y]);
        return { A * k * std::pow(r, k - 1) - B * k * std::pow(r, -k - 1), // strain_rr
                 A * std::pow(r, k - 1) + B * std::pow(r, -k - 1),         // strain_ff
                 T{0} };                                                   // strain_rf
    };

    const auto Expected_Polar_Stress = [&Expected_Polar_Strain, &elastic](const std::array<T, 2>& point) -> std::array<T, 3> {
        static constexpr std::array<T, 2> Zero_Angle_Point = {T{1}, T{0}};
        return calc_stress<T>(elastic.hooke(Zero_Angle_Point), Expected_Polar_Strain(point));
    };

    "displacement"_test = [&mesh, &solution, &Expected_Displacement_R] {
        const auto Expected_Displacement = [&Expected_Displacement_R](const std::array<T, 2>& point) {
            return polar_to_cartesian(Expected_Displacement_R(point), point);
        };
        static constexpr T Epsilon = 8e-4;
        const T error = norm_error(solution.displacement(), mesh->container(), Expected_Displacement);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "strain"_test = [&mesh, &solution, &Expected_Polar_Strain] {
        const auto Expected_Strain = [&Expected_Polar_Strain](const std::array<T, 2>& point) {
            return polar_to_cartesian(Expected_Polar_Strain(point), point);
        };
        static constexpr T Epsilon = 7.3e-2;
        const T error = norm_error(solution.strain(), mesh->container(), Expected_Strain);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress"_test = [&mesh, &solution, &Expected_Polar_Stress] {
        const auto Expected_Stress = [&Expected_Polar_Stress](const std::array<T, 2>& point) {
            return polar_to_cartesian(Expected_Polar_Stress(point), point);
        };
        static constexpr T Epsilon = 4.5e-2;
        const T error = norm_error(solution.stress(), mesh->container(), Expected_Stress);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
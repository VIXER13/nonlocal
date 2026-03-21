#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

#include <embedded_files/composite_ring_su2.h>

namespace {

using T = double;
using I = int64_t;
using namespace boost::ut;
using namespace nonlocal;
using namespace unit_tests;
using namespace mesh;
using namespace solver_2d::mechanical;

constexpr T Expected_Error = T{0};
    
const suite<"isotropic_lame_composite_ring"> _ = [] {
    std::stringstream stream{composite_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T, I>>(stream, mesh_format::SU2);
    const raw_mechanical_parameters<T> parameters = { 
        {"Inner_Material", { .physical = { .elastic = isotropic_elastic_parameters<T>{ .young_modulus = 350., .poissons_ratio = 0.25 } } }},
        {"Outer_Material", { .physical = { .elastic = isotropic_elastic_parameters<T>{ .young_modulus = 150., .poissons_ratio = 0.2  } } }}
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
    const auto solution = equilibrium_equation<I>(mesh, parameters, boundaries_conditions);

    "displacement_x"_test = [&mesh, &solution] {
        static constexpr auto Expected_Displacement_X = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 1.1e-3;
        const T error = norm_error(solution.displacement()[X], mesh->container(), Expected_Displacement_X);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "displacement_y"_test = [&mesh, &solution] {
        static constexpr auto Expected_Displacement_Y = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 1.1e-3;
        const T error = norm_error(solution.displacement()[Y], mesh->container(), Expected_Displacement_Y);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "strain_xx"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_XX = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.1e-2;
        const T error = norm_error(solution.strain()[XX], mesh->container(), Expected_Strain_XX);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "strain_yy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_YY = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.1e-2;
        const T error = norm_error(solution.strain()[YY], mesh->container(), Expected_Strain_YY);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "strain_xy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_XY = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.2e-2;
        const T error = norm_error(solution.strain()[XY], mesh->container(), Expected_Strain_XY);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_xx"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_XX = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.3e-2;
        const T error = norm_error(solution.stress()[XX], mesh->container(), Expected_Stress_XX);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_yy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_YY = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.3e-2;
        const T error = norm_error(solution.stress()[YY], mesh->container(), Expected_Stress_YY);
        //expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_xy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_XY = [](const std::array<T, 2>& point) {
            return T{0};
        };
        static constexpr T Epsilon = 3.3e-2;
        const T error = norm_error(solution.stress()[XY], mesh->container(), Expected_Stress_XY);
        //expect(approx(error, Expected_Error, Epsilon));
    };
};

}
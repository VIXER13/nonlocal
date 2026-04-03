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
            const auto& [x, y] = point;
            const T hypot2 = x * x + y * y;
            if (hypot2 < T{ 0.5625 }) {
                return x * (T{ 2.52841 } + hypot2) / (T{ 47303 } * hypot2);
            }
            else {
                return x * (T{ 9.75721 } - hypot2) / (T{ 140715 } * hypot2);
            }
        };
        static constexpr T Epsilon = 1.6e-3;
        const T error = norm_error(solution.displacement()[X], mesh->container(), Expected_Displacement_X);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "displacement_y"_test = [&mesh, &solution] {
        static constexpr auto Expected_Displacement_X = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T hypot2 = x * x + y * y;
            if (hypot2 < T{ 0.5625 }) {
                return y * (T{ 2.52841 } + hypot2) / (T{ 47303 } * hypot2);
            }
            else {
                return y * (T{ 9.75721 } - hypot2) / (T{ 140715 } * hypot2);
            }
            };
        static constexpr T Epsilon = 1.6e-3;
        const T error = norm_error(solution.displacement()[Y], mesh->container(), Expected_Displacement_X);
        expect(approx(error, Expected_Error, Epsilon));
        };

    "strain_xx"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_XX = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return (T{ 1 } + T{ 2.52841 } * (y2 - x2) / power<2>(hypot2)) / T{ 47303 };
            }
            else {
                return (T{ -1 } + T{ 9.7572 } * (y2 - x2) / power<2>(hypot2)) / T{ 140715 };
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.strain()[XX], mesh->container(), Expected_Strain_XX);
        expect(approx(error, Expected_Error, Epsilon));
        };

    "strain_yy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_YY = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return (T{ 1 } - T{ 2.52841 } * (y2 - x2) / power<2>(hypot2)) / T{ 47303 };
            }
            else {
                return (T{ -1 } - T{ 9.7572 } * (y2 - x2) / power<2>(hypot2)) / T{ 140715 };
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.strain()[YY], mesh->container(), Expected_Strain_YY);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "strain_xy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Strain_XY = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return -x * y / (T{ 9354.31 } * power<2>(hypot2));
            }
            else {
                return -x * y / (T{ 7210.83 } * power<2>(hypot2));
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.strain()[XY], mesh->container(), Expected_Strain_XY);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_xx"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_XX = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return ((power<2>(x2) + power<2>(y2)) / T{ 101.364 } + (y2 - x2) / T{ 66.8165 } + (x2 * y2) / T{ 50.6818 }) 
                    / (power<2>(hypot2));
            }
            else {
                return ((-power<2>(x2) - power<2>(y2)) / T{ 750.481 } + (y2 - x2) / T{ 115.373 } - (x2 * y2) / T{ 375.24 })
                    / (power<2>(hypot2));
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.stress()[XX], mesh->container(), Expected_Stress_XX);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_yy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_YY = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T x2 = x * x;
            const T y2 = y * y;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return ((power<2>(x2) + power<2>(y2)) / T{ 101.364 } - (y2 - x2) / T{ 66.8165 } + (x2 * y2) / T{ 50.6818 })
                    / (power<2>(hypot2));
            }
            else {
                return ((-power<2>(x2) - power<2>(y2)) / T{ 750.481 } - (y2 - x2) / T{ 115.373 } - (x2 * y2) / T{ 375.24 })
                    / (power<2>(hypot2));
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.stress()[YY], mesh->container(), Expected_Stress_YY);
        expect(approx(error, Expected_Error, Epsilon));
    };

    "stress_xy"_test = [&mesh, &solution] {
        static constexpr auto Expected_Stress_XY = [](const std::array<T, 2>& point) {
            const auto& [x, y] = point;
            const T hypot2 = x * x + y * y;
            using metamath::functions::power;
            if (hypot2 < T{ 0.5625 }) {
                return -x * y / (T{ 33.4082 } * power<2>(hypot2));
            }
            else {
                return -x * y / (T{ 57.6866 } * power<2>(hypot2));
            }
            };
        static constexpr T Epsilon = 1.e-1;
        const T error = norm_error(solution.stress()[XY], mesh->container(), Expected_Stress_XY);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
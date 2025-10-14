#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/mechanical/equilibrium_equation_2d.hpp>

#include <boost/ut.hpp>

#include <embedded_files/plate_10x1_h0_125_su2.h>

namespace {

using T = double;
using I = int64_t;
using namespace boost::ut;
using namespace nonlocal;
using namespace mesh;
using namespace solver_2d::mechanical;

std::vector<T> get_values_on_center_line(const mesh_container_2d<T, I>& mesh, const std::vector<T>& solution) {
    std::vector<T> values;
    for(const size_t node : mesh.nodes())
        if (std::abs(mesh.node_coord(node)[X]) < std::numeric_limits<T>::epsilon())
            values.push_back(solution[node]);
    return values;
}

const suite<"saint_venant"> _ = [] {
    std::stringstream stream{plate_10x1_h0_125_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T, I>>(stream, mesh_format::SU2);
    const mechanical_parameters_2d<T> parameters = { 
        .materials = {
            {"DEFAULT", { .physical = { .youngs_modulus = {1, 1}, .poissons_ratio = {0.3, 0.3} } }}
        }
    };
    mechanical_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Vertical"] = {
        std::make_unique<displacement_2d<T>>(T{0}),
        std::make_unique<pressure_2d<T>>(T{0})
    };
    boundaries_conditions["Horizontal"] = {
        std::make_unique<pressure_2d<T>>(T{0}),
        std::make_unique<displacement_2d<T>>(T{0})
    };
    boundaries_conditions["Left"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return -4 * std::abs(point[Y]); }),
        std::make_unique<pressure_2d<T>>(T{0})
    };
    boundaries_conditions["Right"] = {
        std::make_unique<pressure_2d<T>>([](const std::array<T, 2>& point) { return  4 * std::abs(point[Y]); }),
        std::make_unique<pressure_2d<T>>(T{0})
    };
    auto solution = equilibrium_equation<I>(mesh, parameters, boundaries_conditions);
    solution.calc_strain_and_stress();

    "displacement_x_on_center_line"_test = [&mesh, &solution] {
        static constexpr T Expected = 0;
        static constexpr T Epsilon = 2e-14;
        for (const T value : get_values_on_center_line(mesh->container(), solution.displacement()[X]))
            expect(approx(value, Expected, Epsilon));
    };

    "stress_xx_on_center_line"_test = [&mesh, &solution] {
        static constexpr T Expected = 1;
        static constexpr T Epsilon = 1.3e-9;
        for (const T value : get_values_on_center_line(mesh->container(), solution.stress()[0]))
            expect(approx(value, Expected, Epsilon));
    };

    "displacement_x_integral"_test = [&mesh, &solution] {
        static constexpr T Expected = 0;
        static constexpr T Epsilon = 7.5e-13;
        const T integral = mesh::utils::integrate(*mesh, solution.displacement()[X]);
            expect(approx(integral, Expected, Epsilon));
    };

    "stress_xx_integral"_test = [&mesh, &solution] {
        static constexpr T Expected = 10;
        static constexpr T Epsilon = 1e-9;
        const T integral = mesh::utils::integrate(*mesh, solution.stress()[0]);
            expect(approx(integral, Expected, Epsilon));
    };
};

}
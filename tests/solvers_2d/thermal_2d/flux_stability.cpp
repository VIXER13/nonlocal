#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>

#include <boost/ut.hpp>

#include <embedded_files/plate_10x1_h0_125_su2.h>

namespace {

using T = double;
using I = int64_t;
using namespace boost::ut;
using namespace nonlocal;
using namespace mesh;
using namespace solver_2d::thermal;

std::vector<T> get_values_on_center_line(const mesh_container_2d<T, I>& mesh, const std::vector<T>& solution) {
    std::vector<T> values;
    for(const size_t node : mesh.nodes())
        if (std::abs(mesh.node_coord(node)[X]) < std::numeric_limits<T>::epsilon())
            values.push_back(solution[node]);
    return values;
}

const suite<"flux_stability"> _ = [] {
    std::stringstream stream{plate_10x1_h0_125_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T, I>>(stream, mesh_format::SU2);
    const parameters_2d<T> parameters = {{ "DEFAULT", { .physical = { .conductivity = T{1} } } }};
    thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Left"]  = std::make_unique<flux_2d<T>>([](const std::array<T, 2>& point) { return -4 * std::abs(point[Y]); });
    boundaries_conditions["Right"] = std::make_unique<flux_2d<T>>([](const std::array<T, 2>& point) { return  4 * std::abs(point[Y]); });
    auto solution = stationary_heat_equation_solver_2d<I>(mesh, parameters, boundaries_conditions, {});
    solution.calc_flux();

    "temperature_on_center_line"_test = [&mesh, &solution] {
        static constexpr T Expected = 0;
        static constexpr T Epsilon = 2e-14;
        for (const T value : get_values_on_center_line(mesh->container(), solution.temperature()))
            expect(approx(value, Expected, Epsilon));
    };

    "flux_x_on_center_line"_test = [&mesh, &solution] {
        static constexpr T Expected = -1;
        static constexpr T Epsilon = 2e-14;
        for (const T value : get_values_on_center_line(mesh->container(), solution.flux()[X]))
            expect(approx(value, Expected, Epsilon));
    };

    "temperature_integral"_test = [&mesh, &solution] {
        static constexpr T Expected = 0;
        static constexpr T Epsilon = 3.3e-14;
        const T integral = mesh::utils::integrate(*mesh, solution.temperature());
        expect(approx(integral, Expected, Epsilon));
    };

    "flux_x_integral"_test = [&mesh, &solution] {
        static constexpr T Expected = -10;
        static constexpr T Epsilon = 2.2e-14;
        const T integral = mesh::utils::integrate(*mesh, solution.flux()[X]);
        expect(approx(integral, Expected, Epsilon));
    };
};

}
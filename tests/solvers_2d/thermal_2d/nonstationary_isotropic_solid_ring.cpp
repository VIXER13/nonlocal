#include <mesh/mesh_2d/mesh_2d.hpp>
#include <mesh/mesh_2d/mesh_2d_utils.hpp>
#include <mesh/mesh_2d/mesh_container_2d_utils.hpp>
#include <solvers/solver_2d/thermal/nonstationary_heat_equation_solver_2d.hpp>
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

constexpr T Expected_Error = T{0};
constexpr T Time_Step = 0.01;
constexpr uintmax_t Steps_Count = 25;

T exact_solution(const T time, const std::array<T, 2>& point) noexcept {
    const auto& [x, y] = point;
    return std::exp(-time) * (x * x + y * y) * std::cos(2 * std::atan2(y, x));
}

T initial_distribution(const std::array<T, 2>& point) noexcept {
    return exact_solution(0, point);
}

std::function<T(const std::array<T, 2>&)> make_right_part(const T time) {
    return [time](const std::array<T, 2>& point) { return -exact_solution(time, point); };
}

thermal_boundaries_conditions_2d<T> make_boundaries(const T time) {
    thermal_boundaries_conditions_2d<T> boundaries_conditions;
    boundaries_conditions["Inner"] = std::make_unique<temperature_2d<T>>([time](const std::array<T, 2>& point) {
        return exact_solution(time, point);
    });
    boundaries_conditions["Outer"] = std::make_unique<temperature_2d<T>>([time](const std::array<T, 2>& point) {
        return exact_solution(time, point);
    });
    return boundaries_conditions;
}

const suite<"nonstationary_isotropic_solid_ring"> _ = [] {
    std::stringstream stream{solid_ring_su2_data};
    const auto mesh = std::make_shared<mesh_2d<T>>(stream, mesh_format::SU2);
    const raw_thermal_parameters<T> parameters = {
        {"DEFAULT", {.physical = {.conductivity = T{1}, .capacity = T{1}, .density = T{1}}}}
    };
    nonstationary_heat_equation_solver_2d<T> solver{mesh};
    solver.compute(parameters, make_boundaries(0), Time_Step, make_right_part(0), initial_distribution);
    for(const size_t step : std::views::iota(1u, Steps_Count + 1)) {
        const T time = step * Time_Step;
        solver.set_boundaries(make_boundaries(time));
        solver.set_inner_flux(make_right_part(time));
        solver.calc_step();
    }
    const auto solution = solver.solution();

    "temperature"_test = [&mesh, &solution] {
        static constexpr auto Expected_Temperature = [](const std::array<T, 2>& point) {
            return exact_solution(Steps_Count * Time_Step, point);
        };
        static constexpr T Epsilon = 1.3e-3;
        const T error = norm_error(solution.temperature(), mesh->container(), Expected_Temperature);
        expect(approx(error, Expected_Error, Epsilon));
    };
};

}
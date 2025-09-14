#include <solvers/solver_1d/thermal/nonstationary_heat_equation_solver_1d.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::mesh;
using namespace nonlocal::solver_1d::thermal;
using namespace metamath::finite_element;

using T = double;
using I = int64_t;

constexpr size_t Steps = 10 + 1;
constexpr T Time_Step = 0.001;

template<class Exact_Solution>
T calculate_max_error(const mesh_1d<T>& mesh, 
                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature, 
                      Exact_Solution&& exact_solution) {
    T max_error = T{0};
    for (const size_t i : std::ranges::iota_view{0u, size_t(temperature.size())}) 
        if (const T error = std::abs(temperature[i] - exact_solution(mesh.node_coord(i))); error > max_error)
            max_error = error;
    return max_error;
}

T solve_thermal_problem_with_relaxation(const std::shared_ptr<mesh_1d<T>>& mesh,
                                        const parameters_1d<T>& parameters,
                                        const std::function<T(const T, const T)>& exact_solution,
                                        const std::function<T(const T, const T)>& right_part) {
    nonstationary_heat_equation_solver_1d<T, I> solver{mesh, parameters, Time_Step};
    solver.initialize_temperature([&solver, &exact_solution](const T x) { return exact_solution(solver.time(), x); });
    thermal_boundaries_conditions_1d<T> boundaries_conditions = {
        std::make_unique<temperature_1d<T>>(exact_solution(solver.time(), 0)),
        std::make_unique<temperature_1d<T>>(exact_solution(solver.time(), mesh->length()))
    };
    solver.compute(boundaries_conditions);
    for(const size_t step : std::ranges::iota_view{1u, Steps}) {
        const T time = solver.time() + Time_Step;
        boundaries_conditions = {
            std::make_unique<temperature_1d<T>>(exact_solution(time, 0)),
            std::make_unique<temperature_1d<T>>(exact_solution(time, mesh->length()))
        };
        solver.calc_step(boundaries_conditions, [time, &right_part](const T x) { return right_part(time, x); });
    }
    return calculate_max_error(*mesh, solver.temperature(), 
                               [time = solver.time(), &exact_solution](const T x) { return exact_solution(time, x); });
}

const suite<"thermal_nonstationary_1d_relaxation"> _ = [] {
    static constexpr T Epsilon = 5e-5;
    static constexpr size_t Order = 1;
    using quadrature = quadrature_1d<T, gauss, Order>;
    using element_1d = element_1d_integrate<T, lagrangian_element_1d, Order>;

    const auto mesh = std::make_shared<mesh_1d<T>>(
        std::make_unique<element_1d>(quadrature{}), 
        std::vector<segment_data<T>>{{ .length = T{1}, .elements = 100 }});

    "exact_solution"_test = [&mesh] {
        static constexpr auto exact_solution = [](const T t, const T x) constexpr noexcept { return x * x * t * t; };
        static constexpr auto right_part = [](const T t, const T x) constexpr noexcept { return -4 + 4 * std::exp(-t) + 4 * t - 2 * t * t + 2 * t * x * x; };
        const parameters_1d<T> parameters = {{ .physical = {.conductivity = T{1}, .capacity = T{1}, .density = T{1}, .relaxation_time = T{1} } }};
        const T error = solve_thermal_problem_with_relaxation(mesh, parameters, exact_solution, right_part);
        expect(lt(error, Epsilon));
    };
};
}
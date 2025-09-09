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
constexpr T Epsilon = 5e-5;

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

const suite<"thermal_nonstationary_1d_relaxation"> _ = [] {
    static constexpr T Initial_Time = 0;
    constexpr size_t Order = 1;
    using quadrature = quadrature_1d<T, gauss, Order>;
    using element_1d = element_1d_integrate<T, lagrangian_element_1d, Order>;

    const auto mesh = std::make_shared<mesh_1d<T>>(
        std::make_unique<element_1d>(quadrature{}), 
        std::vector<segment_data<T>>{{ .length = T{1}, .elements = 100 }});

    "exact_solution"_test = [&mesh] {
        static constexpr auto exact_solution = [](const T t, const T x) constexpr noexcept { return x * x * t * t; };
        const parameters_1d<T> parameters = {{
            .physical = {
                .conductivity = T{1}, 
                .capacity = T{1},
                .density = T{1},
                .relaxation_time = T{1}
            } 
        }};

        nonstationary_heat_equation_solver_1d<T, I> solver{mesh, parameters, Time_Step};
        solver.initialize_temperature([&solver](const T x) { return exact_solution(solver.time(), x); });
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
            const auto right_part = [t = time](const T x) { return -4 + 4 * std::exp(-t) + 4 * t - 2 * t * t + 2 * t * x * x; };
            solver.calc_step(boundaries_conditions, right_part);
        }

        const T error = calculate_max_error(*mesh, solver.temperature(), 
                                            [time = solver.time()](const T x) { return exact_solution(time, x); });
        expect(lt(error, Epsilon));
    };
};
}
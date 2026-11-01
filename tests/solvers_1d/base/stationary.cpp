#include <solvers/solver_1d/base/stationary_solver_1d.hpp>
#include <tests/solvers_1d/utils/solution_checks.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using I = int64_t;
using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using namespace nonlocal;
using namespace nonlocal::solver_1d;
using namespace nonlocal::unit_tests;
using namespace boost::ut;

static constexpr T Epsilon = T{ 1e-12 };

const suite<"stationary_solver_1d"> _ = [] {
    "neumann_with_volume_load"_test = [] {
        const auto mesh = make_small_mesh<T>();
        // Pass only the model and scalar coefficient to the common solver, without a physical adapter.
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const boundaries_data_1d<T> boundaries{ boundary_data_1d<T>{ .load = T{ -2 } }, boundary_data_1d<T>{ .load = T{ -2 } } };
        const stationary_solver_parameters_1d<T> solver_parameters{
            .right_part = [](const T) constexpr noexcept { return T{ 2 }; }, .integral_value = T{ 8 } / T{ 3 }
        };
        // Volume load 4 balances boundary loads -2 - 2; integral(1 + x - x^2/2) = 8/3.
        // Its piecewise-linear interpolant has integral 5/2, so the discrete constraint shifts all nodes by 1/12.
        const auto settings = init_boundary_settings(problem_settings{ .theories = { theory_t::LOCAL } }, boundaries, true);
        const vector solution = stationary_solver_1d<T, I>(
            mesh, parameters, boundaries, settings, solver_parameters, std::identity{});
        static constexpr auto Expected = std::array{ T{ 13 } / T{ 12 }, T{ 19 } / T{ 12 }, T{ 13 } / T{ 12 }, T{ 0 } };
        check_values(solution, Expected);
    };

    "nonlinear_boundary_iteration"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const boundaries_data_1d<T> boundaries{
            boundary_data_1d<T>{ .prescribed_value = T{ 0 } },
            boundary_data_1d<T>{ .load = T{ 2 },
                .nonlinear = nonlinear_boundary_operator_1d<T>{
                    .value = [](const T u) constexpr noexcept { return -u * u * u; },
                    .derivative = [](const T u) constexpr noexcept { return T{ -3 } * u * u; }
                } }
        };
        stationary_solver_parameters_1d<T> solver_parameters{
            .initial_distribution = [](const T x) constexpr noexcept { return x; }, .tolerance = Epsilon
        };
        const auto settings = init_boundary_settings(problem_settings{ .theories = { theory_t::LOCAL } }, boundaries, true);
        // u(0) = 0, 2*u'(2) = 2 - u(2)^3: u = x/2, since u(2) + u(2)^3 = 2.
        static constexpr auto Expected = std::array{ T{ 0 }, T{ 0.5 }, T{ 1 } };
        check_values(stationary_solver_1d<T, I>(mesh, parameters, boundaries, settings,
                         solver_parameters, std::identity{}), Expected);
        solver_parameters.max_iterations = 1;
        // From u(2) = 2, one Newton step gives 2 - (2 + 8 - 2)/(1 + 12) = 18/13.
        static constexpr auto First_Iteration = std::array{ T{ 0 }, T{ 9 } / T{ 13 }, T{ 18 } / T{ 13 } };
        check_values(stationary_solver_1d<T, I>(mesh, parameters, boundaries, settings,
                         solver_parameters, std::identity{}), First_Iteration);
    };

    "symmetric_indefinite_boundary_jacobian"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const boundaries_data_1d<T> boundaries{
            boundary_data_1d<T>{ .prescribed_value = T{ 0 } },
            boundary_data_1d<T>{ .load = T{ -1 },
                .nonlinear = nonlinear_boundary_operator_1d<T>{
                    .value = [](const T u) constexpr noexcept { return T{ 2 } * u; },
                    .derivative = [](const T) constexpr noexcept { return T{ 2 }; }
                } }
        };
        // 2*u'(2) = -1 + 2*u(2) gives u = x/2. The free-node Jacobian [4 -2; -2 0] is indefinite.
        static constexpr auto Expected = std::array{ T{ 0 }, T{ 0.5 }, T{ 1 } };
        const auto settings = init_boundary_settings(problem_settings{ .theories = { theory_t::LOCAL } }, boundaries, true);
        check_values(stationary_solver_1d<T, I>(mesh, parameters, boundaries, settings, {}, std::identity{}), Expected);
    };
};

} // namespace

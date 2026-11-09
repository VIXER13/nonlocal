#include <solvers/solver_1d/base/boundary_condition_nonlinear_1d.hpp>
#include <solvers/solver_1d/base/stationary_solver_1d.hpp>
#include <solvers/solver_1d/thermal/thermal_boundary_conditions_1d.hpp>
#include <tests/solvers_1d/utils/solution_checks.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using I = int64_t;
using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using sparse_matrix = Eigen::SparseMatrix<T, Eigen::RowMajor, I>;
using namespace nonlocal;
using namespace nonlocal::solver_1d;
using namespace nonlocal::unit_tests;
using namespace boost::ut;

static constexpr T Epsilon = T{ 1e-12 };

const suite<"stationary_solver_1d"> _ = [] {
    "neumann_with_volume_load"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const thermal::thermal_boundaries_conditions_1d<T> boundaries{
            std::make_unique<thermal::flux_1d<T>>(T{ -2 }), std::make_unique<thermal::flux_1d<T>>(T{ -2 })
        };
        const stationary_solver_parameters_1d<T> solver_parameters{
            .right_part = [](const T) constexpr noexcept { return T{ 2 }; }, .integral_value = T{ 8 } / T{ 3 }
        };
        // Volume load 4 balances boundary loads -2 - 2; integral(1 + x - x^2/2) = 8/3.
        // Its piecewise-linear interpolant has integral 5/2, so the discrete constraint shifts all nodes by 1/12.
        const problem_settings settings{ .theories = { theory_t::LOCAL }, .is_neumann = true };
        const vector solution = stationary_solver_1d<T, I>(
            mesh, parameters, boundaries, settings, solver_parameters, std::identity{});
        static constexpr auto Expected = std::array{ T{ 13 } / T{ 12 }, T{ 19 } / T{ 12 }, T{ 13 } / T{ 12 }, T{ 0 } };
        check_values(solution, Expected);
    };

    "nonlinear_boundary_iteration"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const thermal::thermal_boundaries_conditions_1d<T> boundaries{
            std::make_unique<thermal::temperature_1d<T>>(T{ 0 }), std::make_unique<thermal::flux_1d<T>>(T{ 2 })
        };
        static constexpr auto Apply_Cubic = [](sparse_matrix& matrix, vector& residual,
                                               const thermal::thermal_boundaries_conditions_1d<T>&,
                                               const vector& previous, const T time_step) {
            const size_t index = previous.size() - 1;
            const T u = previous[index];
            boundary_condition_nonlinear_1d<true>(matrix, residual, index, u, -u * u * u, T{ -3 } * u * u, time_step);
        };
        stationary_solver_parameters_1d<T> solver_parameters{
            .initial_distribution = [](const T x) constexpr noexcept { return x; }, .tolerance = Epsilon
        };
        const problem_settings settings{
            .theories = { theory_t::LOCAL }, .is_nonlinear_boundary = true, .is_first_kind = { true, false }
        };
        // u(0) = 0, 2*u'(2) = 2 - u(2)^3: u = x/2, since u(2) + u(2)^3 = 2.
        static constexpr auto Expected = std::array{ T{ 0 }, T{ 0.5 }, T{ 1 } };
        check_values(stationary_solver_1d<T, I, physics_t::THERMAL>(
                         mesh, parameters, boundaries, settings, solver_parameters, std::identity{}, nullptr, Apply_Cubic),
                     Expected);
        solver_parameters.max_iterations = 1;
        // From u(2) = 2, one Newton step gives 2 - (2 + 8 - 2)/(1 + 12) = 18/13.
        static constexpr auto First_Iteration = std::array{ T{ 0 }, T{ 9 } / T{ 13 }, T{ 18 } / T{ 13 } };
        check_values(stationary_solver_1d<T, I, physics_t::THERMAL>(
                         mesh, parameters, boundaries, settings, solver_parameters, std::identity{}, nullptr, Apply_Cubic),
                     First_Iteration);
    };

    "symmetric_indefinite_boundary_jacobian"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const std::array parameters{ std::pair{ model_parameters<1, T>{}, coefficient_t<T, 1>{ T{ 2 } } } };
        const thermal::thermal_boundaries_conditions_1d<T> boundaries{
            std::make_unique<thermal::temperature_1d<T>>(T{ 0 }), std::make_unique<thermal::flux_1d<T>>(T{ -1 })
        };
        static constexpr auto Apply_Boundary = [](sparse_matrix& matrix, vector& residual,
                                                  const thermal::thermal_boundaries_conditions_1d<T>&,
                                                  const vector& previous, const T time_step) {
            const size_t index = previous.size() - 1;
            const T u = previous[index];
            boundary_condition_nonlinear_1d<true>(matrix, residual, index, u, T{ 2 } * u, T{ 2 }, time_step);
        };
        // 2*u'(2) = -1 + 2*u(2) gives u = x/2. The free-node Jacobian [4 -2; -2 0] is indefinite.
        static constexpr auto Expected = std::array{ T{ 0 }, T{ 0.5 }, T{ 1 } };
        const problem_settings settings{
            .theories = { theory_t::LOCAL }, .is_nonlinear_boundary = true, .is_first_kind = { true, false }
        };
        check_values(stationary_solver_1d<T, I, physics_t::THERMAL>(
                         mesh, parameters, boundaries, settings, {}, std::identity{}, nullptr, Apply_Boundary),
                     Expected);
    };
};

} // namespace

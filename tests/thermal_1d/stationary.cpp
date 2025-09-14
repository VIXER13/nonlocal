#include "utils.hpp"

#include <metamath/metamath.hpp>
#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_1d/thermal/stationary_heat_equation_solver_1d.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::mesh;
using namespace nonlocal::unit_tests;
using namespace nonlocal::solver_1d::thermal;
using namespace metamath::constants;

template<class T>
using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, std::size_t(1)>;
template<class T>
using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, std::size_t(1)>;

template<std::floating_point T>
void check_solution(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const heat_equation_solution_1d<T>& solution,
                    std::function<T(T)> ref, T eps = T{1e-8}) {
    auto& sol = solution.temperature();
    T err = 0.0, ref_norm = 0.0;
    for (std::size_t k = 0; k < sol.size(); ++k) {
        const T ref_val = ref(mesh->node_coord(k));
        err += std::pow(sol[k] - ref_val, 2);
        ref_norm += std::pow(ref_val, 2);
        //expect(lt(err, eps));
    }
    expect(lt(err / ref_norm, eps));
}

const suite<"thermal_stationary_1d"> _ = [] {
    using T = double;
    using I = int64_t;

    "temperature_temperature"_test = [] {
        // Temperature boundary condition
        // d/dx(k dT/dx) + qv = 0
        // qv(x) = -4 * er * sigma * (-3(x - 2))^(-7/3)
        // T|x=0 = 6^(-1/3), T|x=1 = 3^(-1/3),
        // Exact solution : T(x) = (-3(x - 2))^(-1/3)
        static constexpr auto Expected_Solution = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
        const parameters_1d<T> parameters = {{ .physical = { .conductivity = T{1} } }};
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique<temperature_1d<T>>(Expected_Solution(0.)),
            std::make_unique<temperature_1d<T>>(Expected_Solution(1.))
        };
        const stationary_equation_parameters_1d<T> additional_parameters {
            .right_part = [](const T x) {  return -4. / (9 * metamath::functions::power<2>(2 - x) * std::cbrt(6. - 3. * x)); },
        };

        T prev_error = std::numeric_limits<T>::max();
        for(const size_t elements : {10, 20, 40}) {
            const auto mesh = std::make_shared<mesh_1d<T>>(
                std::make_unique<element_1d<T>>(quadrature<T>{}),
                std::vector<segment_data<T>>{{ .length = T{1}, .elements = elements }});
            const auto solution = stationary_heat_equation_solver_1d<T, I>(mesh, parameters, boundaries_conditions, additional_parameters);
            const auto expected_solution_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Solution);
            const T error = max_error(solution.temperature(), expected_solution_discrete) / max_norm(expected_solution_discrete);
            expect(lt(error, prev_error));
            prev_error = error;
        }
    };

    "const_temperature_radiation"_test = [] {
        // Radiation boundary condition
        // d/dx(lambda * dT/dx) + qv = 0
        // alpha = er * sigma / lambda
        // qv(x) = -4 * lambda * alpha * alpha * (-3 * alpha * (x - 1.01))^(-7/3)
        // T|x=0 = T_ref(0), q*n|x=L = er * sigma * T^4,
        // Exact solution : T_ref(x) = (-3 * alpha * (x - 1.01))^(-1/3)
        constexpr T sigma = Stefan_Boltzmann_Constant<T>;
        constexpr T lambda = 10.;
        constexpr T emissivity = 0.7;
        constexpr T alpha = emissivity * sigma / lambda;
        const auto ref_sol = [&](T x) constexpr noexcept -> T {
            return std::pow(-3. * alpha * (x - 1.01), -1./3.);
        };
        constexpr T left_temperature = ref_sol(0.0);
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique<temperature_1d<T>>(left_temperature),
            std::make_unique<radiation_1d<T>>(emissivity)
        };
        const std::vector<mesh::segment_data<T>> segments({{ .length = T(1.0), .search_radius = T(0.0), .elements = I(1000) }});
        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .physical = { .conductivity = lambda } };
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(std::make_unique<element_1d<T>>(quadrature<T>()), segments);
        const stationary_equation_parameters_1d<T> additional_parameters {
            .right_part = [&](const T x) constexpr noexcept {  return -4. * lambda * alpha * alpha * std::pow(-3. * alpha * (x - 1.01), -7./3.); },
            .initial_distribution = [ref_sol](const T x) constexpr noexcept { return ref_sol(1.0); },
            .tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-8,
            .max_iterations = 10,
            .energy = T{0}
        };
        const auto num_sol = stationary_heat_equation_solver_1d<T, I>(mesh, parameters, boundaries_conditions, additional_parameters);
        check_solution<T>(mesh, num_sol, ref_sol, 1e-6);
    };

    "radiation_const_temperature"_test = [] {
        // Radiation boundary condition
        // d/dx(lambda * dT/dx) + qv = 0
        // alpha = er * sigma / lambda
        // qv(x) = -4 * lambda * alpha * alpha * (3 * alpha * (x - 1.01))^(-7/3)
        // q*n|x=0 = = er * sigma * T^4, T|x=L = T_ref(L),
        // Exact solution : T_ref(x) = (3 * alpha * (x - 1.01))^(-1/3)
        constexpr T sigma = Stefan_Boltzmann_Constant<T>;
        constexpr T lambda = 10.;
        constexpr T emissivity = 0.7;
        constexpr T alpha = emissivity * sigma / lambda;
        const auto ref_sol = [&](T x) constexpr noexcept -> T {
            return std::pow( T(3.0) * alpha * (x + T(1.01)), T(-1./3.) );
        };
        constexpr T right_temperature = ref_sol(1.0);
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique<radiation_1d<T>>(emissivity),
            std::make_unique<temperature_1d<T>>(right_temperature)
        };
        const std::vector<mesh::segment_data<T>> segments({{ .length = T(1.0), .search_radius = T(0.0), .elements = I(1000) }});
        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .physical = { .conductivity = lambda } };
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(std::make_unique<element_1d<T>>(quadrature<T>()), segments);
        const stationary_equation_parameters_1d<T> additional_parameters {
            .right_part = [&](const T x) constexpr noexcept {  return -4. * lambda * alpha * alpha * std::pow(3. * alpha * (x + 1.01), -7./3.); },
            .initial_distribution = [ref_sol](const T x) constexpr noexcept { return ref_sol(0.0) + 1.; },
            .tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-8,
            .max_iterations = 10,
            .energy = T{0}
        };
        const auto num_sol = stationary_heat_equation_solver_1d<T, I>(mesh, parameters, boundaries_conditions, additional_parameters);
        check_solution<T>(mesh, num_sol, ref_sol, 1e-6);
    };
};
    
}

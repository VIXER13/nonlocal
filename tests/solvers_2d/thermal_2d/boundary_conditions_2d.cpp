#include <solvers/solver_2d/thermal/stationary_heat_equation_solver_2d.hpp>

#include <boost/ut.hpp>
#include <nlohmann/json.hpp>

#include <embedded_files/square_1x1_16el_uniform_mesh_su2.h>

namespace {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::solver_2d::thermal;
using namespace metamath::finite_element;
using namespace metamath::constants;

template<std::floating_point T, std::signed_integral I>
void check_solution(const std::shared_ptr<mesh::mesh_2d<T, I>>& mesh, const heat_equation_solution_2d<T, I>& solution,
                    std::function<T(const std::array<T, 2>&)> ref, T eps = T{1e-8}) {
    auto& sol = solution.temperature();
    const auto& container = mesh->container();
    T err = 0.0, ref_norm = 0.0;
    for (std::size_t k = 0; k < sol.size(); ++k) {
        const T ref_val = ref(container.node_coord(k));
        err += std::pow(sol[k] - ref_val, 2);
        ref_norm += std::pow(ref_val, 2);
    }
    expect(lt(err / ref_norm, eps));
}

template<std::floating_point T, std::signed_integral I>
std::tuple<std::shared_ptr<mesh::mesh_2d<T, I>>, parameters_2d<T>,
           stationary_equation_parameters_2d<T>, 
           std::function<T(const std::array<T, 2>&)>> test_task_1(double mult = 1.0) {
    // div(k grad T) + qv = kx d^2T/dx^2 + ky d^2T/dy^2 + qv = 0
    // Exact solution : T(x, y) = f(x) * g(y) 
    // Let f be : df/dx = alpha * f^4
    // f(x) = (-3 * alpha * (x - 1.01))^(-1/3)
    // g(y) = A * exp( -(y - m)^2 / d )
    // qv(x, y) = -(kx * f'' * g +  ky * f * g'')
    // Temperature BC
    // T(0, y) = f(0) * g(y), T(Lx, y) = f(Lx) * g(y)
    // T(x, 0) = f(x) * g(0), T(x, Ly) = f(x)  * g(Ly)
    // Radiation BC
    // q*n|x=Lx = er * sigma * T^4
    constexpr T sigma = Stefan_Boltzmann_Constant<T>;
    constexpr T lambda = 10.;
    constexpr T emissivity = 0.7;
    const T alpha = mult * emissivity * sigma / lambda;
    const T shift = mult * 1.01;
    const auto f = [alpha, shift](T x) constexpr noexcept -> T {
        return std::pow(-3. * alpha * (x - shift), -1./3.);
    };
    constexpr T A = 1., m = 0.5, d = 0.1;
    const auto g = [](T y) constexpr noexcept -> T {
        return A * std::exp( -(y - m) * (y - m) / d );
    };
    const auto ref_sol = [f, g](const std::array<T, 2>& x) constexpr noexcept -> T {
        return f(x[0]) * g(x[1]);
    };
    // Computational mesh
    std::stringstream stream{square_1x1_16el_uniform_mesh_su2_data};
    const auto mesh = std::make_shared<nonlocal::mesh::mesh_2d<double, int64_t>>(stream, mesh::mesh_format::SU2);
    // Material and model parameters
    parameters_2d<T> parameters;
    parameters["Layer1"] = {
        .physical = {
            .conductivity = anisotropic_conductivity_t<T> {lambda, 0.0, 0.0, lambda},
            .capacity = T{1},
            .density = T{1},
            .relaxation_time = T{0}
        }
    };
    // Auxilary data
    const auto right_part = [ref_sol, alpha, f](const std::array<T, 2>& x) constexpr noexcept -> T { 
        const T a1 = 4 * alpha * alpha * lambda;
        constexpr T a2 = 2 * lambda / d;
        return -ref_sol(x) * (a1 * metamath::functions::power<6>(f(x[0])) + a2 * (2. / d * metamath::functions::power<2>(x[1] - m) - 1.)); 
    };
    constexpr auto initial_distribution = [&](const std::array<T, 2>& x) constexpr noexcept -> T { 
        return 2.; 
    };
    const stationary_equation_parameters_2d<T> auxiliary_data {
        .right_part = right_part,
        .initial_distribution = initial_distribution,
        .tolerance = std::is_same_v<T, float> ? 1e-5 : 1e-10,
        .max_iterations = 10,
        .energy = T{0.0}
    };
    return std::make_tuple(mesh, parameters, auxiliary_data, ref_sol);
}

const suite<"thermal_stationary_boundary_conditions_2d"> _ = [] {
    using T = double;
    using I = int64_t;

    "all_boundaries_temperature"_test = [] {
        const auto task = test_task_1<T, I>();
        auto& mesh = std::get<0>(task);
        auto& parameters = std::get<1>(task);
        auto& auxiliary_data = std::get<2>(task);
        auto& ref_sol = std::get<3>(task);
        // Boundaries conditions
        const auto left_temperature   = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({0.0, x[1]}); };
        const auto right_temperature  = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({1.0, x[1]}); };
        const auto top_temperature    = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 1.0}); };
        const auto bottom_temperature = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 0.0}); };
        thermal_boundaries_conditions_2d<T> boundaries_conditions;
        boundaries_conditions["Left"  ] = std::make_unique<temperature_2d<T>>(left_temperature);
        boundaries_conditions["Right" ] = std::make_unique<temperature_2d<T>>(right_temperature);
        boundaries_conditions["Top"   ] = std::make_unique<temperature_2d<T>>(top_temperature);
        boundaries_conditions["Bottom"] = std::make_unique<temperature_2d<T>>(bottom_temperature);
        // 2D stationary solution
        const auto num_sol = stationary_heat_equation_solver_2d<I, T, I>( 
            mesh, parameters, boundaries_conditions, auxiliary_data
        );
        check_solution<T, I>(mesh, num_sol, ref_sol, 1e-3);
    };

    // "radiation_on_right_side"_test = [] {
    //     const auto [mesh, parameters, auxiliary_data, ref_sol] = test_task_1<T, I>(-1.0);
    //     constexpr T emissivity = 0.7;
    //     // Boundaries conditions
    //     const auto left_temperature   = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({0.0, x[1]}); };
    //     const auto top_temperature    = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 1.0}); };
    //     const auto bottom_temperature = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 0.0}); };
    //     thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
    //     boundaries_conditions["Left"  ] = std::make_unique<thermal::temperature_2d<T>>(left_temperature);
    //     boundaries_conditions["Right" ] = std::make_unique<thermal::radiation_2d<T>>(emissivity);
    //     boundaries_conditions["Top"   ] = std::make_unique<thermal::temperature_2d<T>>(top_temperature);
    //     boundaries_conditions["Bottom"] = std::make_unique<thermal::temperature_2d<T>>(bottom_temperature);
    //     // 2D stationary solution
    //     const auto num_sol = nonlocal::thermal::stationary_heat_equation_solver_2d<I, T, I>( 
    //         mesh, parameters, boundaries_conditions, auxiliary_data
    //     );
    //     check_solution<T, I>(mesh, num_sol, ref_sol, 1e-3);
    // };

    // "radiation_on_left_side"_test = [] {
    //     const auto [mesh, parameters, auxiliary_data, ref_sol] = test_task_1<T, I>();
    //     constexpr T emissivity = 0.7;
    //     // Boundaries conditions
    //     const auto right_temperature  = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({1.0, x[1]}); };
    //     const auto top_temperature    = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 1.0}); };
    //     const auto bottom_temperature = [&](const std::array<T, 2>& x) constexpr noexcept -> T { return ref_sol({x[0], 0.0}); };
    //     thermal::thermal_boundaries_conditions_2d<T> boundaries_conditions;
    //     boundaries_conditions["Left"  ] = std::make_unique<thermal::radiation_2d<T>>(emissivity);
    //     boundaries_conditions["Right" ] = std::make_unique<thermal::temperature_2d<T>>(right_temperature);
    //     boundaries_conditions["Top"   ] = std::make_unique<thermal::temperature_2d<T>>(top_temperature);
    //     boundaries_conditions["Bottom"] = std::make_unique<thermal::temperature_2d<T>>(bottom_temperature);
    //     // 2D stationary solution
    //     const auto num_sol = nonlocal::thermal::stationary_heat_equation_solver_2d<I, T, I>( 
    //         mesh, parameters, boundaries_conditions, auxiliary_data
    //     );
    //     check_solution<T, I>(mesh, num_sol, ref_sol, 1e-3);
    // };

};
    
}

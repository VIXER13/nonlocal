#include "metamath.hpp"
#include "logger.hpp"
#include "thermal/stationary_heat_equation_solver_1d.hpp"
#include "thermal/nonstationary_heat_equation_solver_1d.hpp"
#include "influence_functions_1d.hpp"

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

template<class T>
using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, std::size_t(1)>;
template<class T>
using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, std::size_t(1)>;

template <std::floating_point T>
struct time_data final {
    T time_step = T{0};       
    T initial_time = T{0};
    uint64_t steps_count = 0; 
};

template<std::floating_point T, std::signed_integral I>
void check_solution(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const heat_equation_solution_1d<T>& solution, T time_layer, I step, 
                    std::function<T(T, T)> ref, T eps = T{1e-8}) {
    auto& sol = solution.temperature();
    for (std::size_t k = 0; k < sol.size(); ++k) {
        expect(lt(std::abs(sol[k] - ref(time_layer, mesh->node_coord(k))), eps * step));
    }
}

template<template<class> class Left_bc_type, template<class> class Right_bc_type, std::floating_point T, std::signed_integral I>
void solve_nonstationary_thermal_1d_problem(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const time_data<T>& time,
                                            const parameters_1d<T>& parameters, const std::function<T(T)>& init_dist, 
                                            const std::function<T(T, T)>& right_part,
                                            const std::function<T(T)>& left_bc, const std::function<T(T)>& right_bc,
                                            std::function<T(T, T)> ref_sol, T eps = T{1e-8}, bool internal_iters_check = false)  {

    nonstationary_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};

    {   // Step initial
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique<Left_bc_type<T>>(left_bc(T(0))),
            std::make_unique<Right_bc_type<T>>(right_bc(T(0)))
        };
        solver.compute(parameters, boundaries_conditions, init_dist);
    }

    for(const I step : std::ranges::iota_view{1u, time.steps_count + 1}) {
        const T time_layer = step * time.time_step; 
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique<Left_bc_type<T>>(left_bc(T(time_layer))),
            std::make_unique<Right_bc_type<T>>(right_bc(T(time_layer)))
        };
        solver.compute(parameters, boundaries_conditions, EMPTY_FUNCTION);  
        const auto rp = [&right_part, &time_layer](const T x) { return right_part(time_layer, x); };
        solver.calc_step(boundaries_conditions, rp);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        if(step == time.steps_count) 
            check_solution<T, I>(mesh, solution, time_layer, step, ref_sol, eps);
        else if(internal_iters_check) 
            check_solution<T, I>(mesh, solution, time_layer, step, ref_sol, eps);
    }
}

const suite<"thermal_nonstationary_boundary_conditions_1d"> _ = [] {

    using T = double;
    using I = int64_t;

    "const_temperature_temperature"_test = [] {
        // Two first kind BCs
        // BC T|x=0 = alpha, T|x=L = beta
        // Exact solution T(x) = alpha * (x - L) / L + beta * x / L
        constexpr T left_temp = T(10.0), right_temp = T(-3.0);
        const auto left_bc  = [&left_temp] (const T t) constexpr noexcept { return left_temp; };
        const auto right_bc = [&right_temp](const T t) constexpr noexcept { return right_temp; };
        const auto ref_sol  = [](T t, T x) {
            return left_temp * (2.0 - x) * 0.5 + right_temp * x * 0.5;
        };

        const std::vector<mesh::segment_data<T>> segments = {{ .length = T(2.0), .search_radius = T(0.0), .elements = I(100) }};

        parameters_1d<T> parameters ({{ .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                       // conductivity, capacity, density, relaxation_time
                                        .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(T(1.0), T(1.0), T(1.0), T(0.0)) }});
        
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(std::make_unique<element_1d<T>>(quadrature<T>()), segments);
        
        constexpr auto init_dist  = [](const T x)            constexpr noexcept { return T(0.0); };
        constexpr auto right_part = [](const T t, const T x) constexpr noexcept { return T(0.0); };
        
        const time_data<T> time ({.time_step = T(0.5), .initial_time = T(0.0), .steps_count = I(10)});

        solve_nonstationary_thermal_1d_problem<temperature_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_bc, right_bc, ref_sol, 1e-3);
    };

    "const_temperature_flux"_test = [] {
        // Zero flux with temperature
        // T|x=0 = alpha, q*n|x=L = 0 OR q*n|x=0 = 0, T|x=L = alpha
        // Exact solution T(x) = alpha
        constexpr T temp = T(10.0), flux = T(0.0);
        const auto temp_bc = [&temp](const T t) constexpr noexcept { return temp; };
        const auto flux_bc = [&flux](const T t) constexpr noexcept { return flux; };
        constexpr auto ref_sol = [](T t, T x) {
            return temp;
        };

        const std::vector<mesh::segment_data<T>> segments({{ .length = T(2.0), .search_radius = T(0.0), .elements = I(100) }});

        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                     // conductivity, capacity, density, relaxation_time
                                     .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(T(1.0), T(1.0), T(1.0), T(0.0)) };
        
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(std::make_unique<element_1d<T>>(quadrature<T>()), segments);

        constexpr auto init_dist  = [](const T x)            constexpr noexcept { return T(0.0); };
        constexpr auto right_part = [](const T t, const T x) constexpr noexcept { return T(0.0); };
        
        const time_data<T> time ({.time_step = T(10.0), .initial_time = T(0.0), .steps_count =  I(10) });

        // T|x=0 = alpha, q*n|x=L = 0
        solve_nonstationary_thermal_1d_problem<temperature_1d, flux_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                              temp_bc, flux_bc, ref_sol, 1e-8);
        // q*n|x=0 = 0, T|x=L = alpha
        solve_nonstationary_thermal_1d_problem<flux_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                              flux_bc, temp_bc, ref_sol, 1e-8);
    };

    "all_nonconst_temperature_and_flux_configurations"_test = [] {
        // General thermal conductivity equation solution
        // T(x,t) = exp(-k * t) Cos(w1 * x) Sin(w2 * x)
        constexpr T k = T(1.0), w1 = T(2.0), w2 = T(3.0);
        constexpr T lambda = T(5.0), cap = T(1000.0), rho = T(1000.0), L = T(2.0);

        constexpr auto left_temp =  [&](const T t) constexpr noexcept { return T(0.0); };
        const auto right_temp = [&](const T t) constexpr noexcept { return std::exp(-k * t) * std::cos(w1 * L) * std::sin(w2 * L); };

        const auto left_flux =  [&](const T t) constexpr noexcept { return -lambda * std::exp(-k * t) * w2; };
        const auto right_flux = [&](const T t) constexpr noexcept { 
            return  lambda * std::exp(-k * t) * (-w1 * std::sin(w1 * L) * std::sin(w2 * L) + w2 * std::cos(w1 * L) * std::cos(w2 * L)); 
        };

        const auto ref_sol = [&](T t, T x) {
            return std::exp(-k * t) * std::cos(w1 * x) * std::sin(w2 * x);
        };

        const std::vector<mesh::segment_data<T>> segments({{ .length = T(L), .search_radius = T(0.0), .elements = I(200) }});

        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                     // conductivity, capacity, density, relaxation_time
                                     .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(lambda, cap, rho, T(0.0)) };
        
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(std::make_unique<element_1d<T>>(quadrature<T>()), segments);
        
        const auto init_dist =  [&](const T x)            constexpr noexcept { return ref_sol(T(0), x); };
        const auto right_part = [&](const T t, const T x) constexpr noexcept { 
            return (-k * cap * rho + w1 * w1* lambda + w2 * w2* lambda) * ref_sol(t, x) + 2 * lambda * w1 * w2 * std::exp(-k * t) * std::sin(w1 * x) * std::cos(w2 * x);
        };
        
        const time_data<T> time ({ .time_step = T(0.005), .initial_time = T(0.0), .steps_count =  I(10) });

        solve_nonstationary_thermal_1d_problem<temperature_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_temp, right_temp, ref_sol, 1e-2, true);
        solve_nonstationary_thermal_1d_problem<temperature_1d, flux_1d,        T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_temp, right_flux, ref_sol, 1e-2, true);
        solve_nonstationary_thermal_1d_problem<flux_1d,        temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_flux, right_temp, ref_sol, 1e-2, true);
    };
};
    
}

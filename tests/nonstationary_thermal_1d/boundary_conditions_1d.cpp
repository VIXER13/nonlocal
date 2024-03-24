#include "nonstationary_tests_1d.hpp"

namespace bc_1d_tests {

using namespace nonstat_1d_tests;
using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::thermal;

template<template<class> class left_bc_type, template<class> class right_bc_type, std::floating_point T, std::signed_integral I>
void solve_nonstationary_thermal_1d_problem(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const config::time_data<T>& time,
                                            const parameters_1d<T>& parameters, std::function<T(T)> init_dist, std::function<T(T, T)> right_part,
                                            std::function<T(T)> LEFT_BC, std::function<T(T)> RIGHT_BC,
                                            std::function<T(T, T)> ref_sol, T eps = epsilon, bool internal_iters_check = false)  {

   nonstationary_heat_equation_solver_1d<T, I> solver{mesh, time.time_step};

    {   // Step initial
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique <left_bc_type<T>>(LEFT_BC(T(0))),
            std::make_unique<right_bc_type<T>>(RIGHT_BC(T(0)))
        };
        solver.compute(parameters, boundaries_conditions, init_dist);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, 0, solution);
    }

    for(const I step : std::ranges::iota_view{1u, time.steps_count + 1}) {
        const T time_layer = (step - 1) * time.time_step; 
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique <left_bc_type<T>>(LEFT_BC(T(time_layer))),
            std::make_unique<right_bc_type<T>>(RIGHT_BC(T(time_layer)))
        };
        solver.compute(parameters, boundaries_conditions, EMPTY_FUNCTION);  
        std::function<T(T)> rp = [&right_part, &time_layer, &time](const T x) { return right_part(time_layer, x); };
        solver.calc_step(boundaries_conditions, rp);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, step, solution);
        if(internal_iters_check) check_solution<T, I>(mesh, solution, time_layer, step, ref_sol, eps);
    }

    {   // Step final
        const T time_layer = time.steps_count * time.time_step; 
        const thermal_boundaries_conditions_1d<T> boundaries_conditions = {
            std::make_unique <left_bc_type<T>>(LEFT_BC(T(time_layer))),
            std::make_unique<right_bc_type<T>>(RIGHT_BC(T(time_layer)))
        };
        solver.compute(parameters, boundaries_conditions, EMPTY_FUNCTION); 
        std::function<T(T)> rp = [&right_part, &time_layer, &time](const T x) { return right_part(time_layer, x); };
        solver.calc_step(boundaries_conditions, rp);
        heat_equation_solution_1d<T> solution{mesh, parameters, solver.temperature()};
        save_and_calc_flux<T, I>(time, time.steps_count + 1, solution);
        check_solution<T, I>(mesh, solution, time_layer, time.steps_count + 1, ref_sol, eps);
    }
}

const suite<"1d_boundary_conditions_temperature_flux"> _ = [] {
    "Test №1"_test = [] {
        // Two first kind BCs
        // BC T|x=0 = alpha, T|x=L = beta
        // Exact solution T(x) = alpha * (x - L) / L + beta * x / L
        constexpr T left_temp = T(10.0), right_temp = T(-3.0);
        std::function<T(T)> left_bc  = [&left_temp] (const T t) constexpr noexcept { return left_temp; };
        std::function<T(T)> right_bc = [&right_temp](const T t) constexpr noexcept { return right_temp; };
        std::function<T(T, T)> ref_sol = [](T t, T x) {
            return left_temp * (2.0 - x) * 0.5 + right_temp * x * 0.5;
        };

        //Same calculation with json file
        // std::filesystem::path path = "PATH_TO_REPOSITORY/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test1.json";
        // solve_nonstationary_thermal_1d_problem_json<T, I>(nonlocal::config::parse_json(path), ref_sol);

        // Set computational domains
        std::vector<mesh::segment_data<T>> segments(I(1));
        segments[0] = mesh::segment_data<T>{ .length = T(2.0), .search_radius = T(0.0), .elements = I(100) };
        // Domain parameters
        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                     // conductivity, capacity, density, relaxation_time
                                     .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(T(1.0), T(1.0), T(1.0), T(0.0)) };
        
        // segments, element_order, quadrature_order 
        const auto mesh = make_mesh_1d(segments, config::order_t::LINEAR, config::order_t::LINEAR);
        
        //right_part, initial_distribution
        std::function<T(T)>     init_dist = [](const T x)            constexpr noexcept { return T(0.0); };
        std::function<T(T, T)> right_part = [](const T t, const T x) constexpr noexcept { return T(0.0); };
        
        // time_step, initial_time, steps_count, save_frequency
        const config::time_data<T> time { T(0.1), T(0.0), I(100), I(1) };

        solve_nonstationary_thermal_1d_problem<temperature_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_bc, right_bc, ref_sol, 1e-8);
    };

    "Test №2"_test = [] {
        // Zero flux with temperature
        // T|x=0 = alpha, q*n|x=L = 0 OR q*n|x=0 = 0, T|x=L = alpha
        // Exact solution T(x) = alpha
        constexpr T temp = T(10.0), flux = T(0.0);
        std::function<T(T)> temp_bc = [&temp](const T t) constexpr noexcept { return temp; };
        std::function<T(T)> flux_bc = [&flux](const T t) constexpr noexcept { return flux; };
        std::function<T(T, T)> ref_sol = [](T t, T x) {
            return temp;
        };

        // Same calculation with json file
        // std::filesystem::path path1 = "PATH_TO_REPOSITORY/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test2_1.json";
        // std::filesystem::path path2 = "PATH_TO_REPOSITORY/nonlocal/tests/nonstationary_thermal_1d/json_files/temperature_flux/test2_2.json";
        // solve_nonstationary_thermal_1d_problem_json<T, I>(nonlocal::config::parse_json(path1), ref_sol);
        // solve_nonstationary_thermal_1d_problem_json<T, I>(nonlocal::config::parse_json(path2), ref_sol);

        // Set computational domains
        std::vector<mesh::segment_data<T>> segments(I(1));
        segments[0] = mesh::segment_data<T>{ .length = T(2.0), .search_radius = T(0.0), .elements = I(100) };
        // Domain parameters
        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                     // conductivity, capacity, density, relaxation_time
                                     .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(T(1.0), T(1.0), T(1.0), T(0.0)) };
        
        // segments, element_order, quadrature_order 
        const auto mesh = make_mesh_1d(segments, config::order_t::LINEAR, config::order_t::LINEAR);
        
        //right_part, initial_distribution
        std::function<T(T)>     init_dist = [](const T x)            constexpr noexcept { return T(0.0); };
        std::function<T(T, T)> right_part = [](const T t, const T x) constexpr noexcept { return T(0.0); };
        
        // time_step, initial_time, steps_count, save_frequency
        const config::time_data<T> time { T(0.1), T(0.0), I(1000), I(1) };

        // T|x=0 = alpha, q*n|x=L = 0
        solve_nonstationary_thermal_1d_problem<temperature_1d, flux_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                              temp_bc, flux_bc, ref_sol, 1e-8);
        // q*n|x=0 = 0, T|x=L = alpha
        solve_nonstationary_thermal_1d_problem<flux_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                              flux_bc, temp_bc, ref_sol, 1e-8);
    };

    "Test №3"_test = [] {
        // General thermal conductivity equation solution
        // T(x,t) = exp(-k * t) Cos(w1 * x) Sin(w2 * x)
        constexpr T k = T(1.0), w1 = T(2.0), w2 = T(3.0);
        constexpr T lambda = T(5), cap = T(1000), rho = T(1000), L = T(2.0);

        std::function<T(T)> left_temp =  [&](const T t) constexpr noexcept { return T(0); };
        std::function<T(T)> right_temp = [&](const T t) constexpr noexcept { return std::exp(-k * t) * std::cos(w1 * L) * std::sin(w2 * L); };

        std::function<T(T)> left_flux =  [&](const T t) constexpr noexcept { return -lambda * std::exp(-k * t) * w2; };
        std::function<T(T)> right_flux = [&](const T t) constexpr noexcept { 
            return  lambda * std::exp(-k * t) * (-w1 * std::sin(w1 * L) * std::sin(w2 * L) + w2 * std::cos(w1 * L) * std::cos(w2 * L)); 
        };

        std::function<T(T, T)> ref_sol = [&](T t, T x) {
            return std::exp(-k * t) * std::cos(w1 * x) * std::sin(w2 * x);
        };

        // Set computational domains
        std::vector<mesh::segment_data<T>> segments(I(1));
        segments[0] = mesh::segment_data<T>{ .length = L, .search_radius = T(0.0), .elements = I(200) };
        // Domain parameters
        parameters_1d<T> parameters(segments.size());
        parameters[0] = { .model = { .influence = influence::polynomial_1d<T, 1, 1>{T(0.0)}, .local_weight = T(1.0) },
                                     // conductivity, capacity, density, relaxation_time
                                     .physical = std::make_shared<parameter_1d<T, coefficients_t::CONSTANTS>>(lambda, cap, rho, T(0.0)) };
        
        // segments, element_order, quadrature_order 
        const auto mesh = make_mesh_1d(segments, config::order_t::LINEAR, config::order_t::LINEAR);
        
        //right_part, initial_distribution
        std::function<T(T)>     init_dist = [&](const T x)            constexpr noexcept { return ref_sol(T(0), x); };
        std::function<T(T, T)> right_part = [&](const T t, const T x) constexpr noexcept { 
            return (-k * cap * rho + w1 * w1* lambda + w2 * w2* lambda) * ref_sol(t, x) + 2 * lambda * w1 * w2 * std::exp(-k * t) * std::sin(w1 * x) * std::cos(w2 * x);
        };
        
        // time_step, initial_time, steps_count, save_frequency
        const config::time_data<T> time { T(1e-4), T(0.0), I(100), I(1) };

        solve_nonstationary_thermal_1d_problem<temperature_1d, temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_temp, right_temp, ref_sol, 1e-3, true);
        solve_nonstationary_thermal_1d_problem<temperature_1d, flux_1d,        T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_temp, right_flux, ref_sol, 1e-3, true);
        solve_nonstationary_thermal_1d_problem<flux_1d,        temperature_1d, T, I>(mesh, time, parameters, init_dist, right_part, 
                                                                                     left_flux, right_temp, ref_sol, 1e-3, true);
    };
};
    
}

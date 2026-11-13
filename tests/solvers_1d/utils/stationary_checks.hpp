#pragma once

#include "solution_checks.hpp"

#include <solvers/solver_1d/base/stationary_solver_1d.hpp>

namespace nonlocal::unit_tests {

// The same examples exercise each physical interface; Solve maps common coefficients and settings to that interface.
template<std::floating_point T, physics_t Physics, template<class> class Prescribed_Value,
         template<class> class Load, template<class> class Robin, class Solve>
void register_stationary_examples(const Solve& solve) {
    using namespace boost::ut;
    using boundaries = nonlocal::solver_1d::boundaries_conditions_1d<T, Physics>;
    using solver_parameters = nonlocal::solver_1d::stationary_solver_parameters_1d<T>;

    "dirichlet_and_volume_load"_test = [solve] {
        const boundaries conditions{
            std::make_unique<Prescribed_Value<T>>(T{ 1 }), std::make_unique<Prescribed_Value<T>>(T{ 1 })
        };
        const solver_parameters parameters{ .right_part = [](const T) constexpr noexcept { return T{ 2 }; } };
        // -2*u'' = 2 on [0, 2], u(0) = u(2) = 1: u = 1 + x - x^2/2.
        static constexpr auto Expected = std::array{ T{ 1 }, T{ 1.5 }, T{ 1 } };
        check_values(solve(make_small_mesh<T>(), T{ 2 }, {}, conditions, parameters), Expected);
    };

    "robin_and_load"_test = [solve] {
        const boundaries conditions{
            std::make_unique<Robin<T>>(T{ 3 }, T{ 1 } / T{ 3 }), std::make_unique<Load<T>>(T{ 2 })
        };
        // u = 1 + x: at x = 0, -2*u' + 3*u = 1; at x = 2, 2*u' = 2.
        static constexpr auto Expected = std::array{ T{ 1 }, T{ 2 }, T{ 3 } };
        check_values(solve(make_small_mesh<T>(), T{ 2 }, {}, conditions, {}), Expected);
    };

    "neumann_integral_value"_test = [solve] {
        const auto mesh = make_small_mesh<T>();
        for (const T integral : { T{ 0 }, T{ 8 } })
            for (const T load : { T{ 0 }, T{ 2 } }) {
                const boundaries conditions{ std::make_unique<Load<T>>(-load), std::make_unique<Load<T>>(load) };
                // 2*u' = load and integral(u) = Q: u(x) = Q/2 + load*(x - 1)/2.
                // This also checks a constant field with nonzero integral and zero loads.
                const auto expected = std::array{ (integral - load) / T{ 2 }, integral / T{ 2 }, (integral + load) / T{ 2 } };
                check_values(solve(mesh, T{ 2 }, {}, conditions, { .integral_value = integral }), expected);
            }
    };

    "incompatible_neumann_loads"_test = [solve] {
        const boundaries conditions{ std::make_unique<Load<T>>(T{ 1 }), std::make_unique<Load<T>>(T{ 0 }) };
        // There is no volume load to balance the total boundary load 1.
        expect(throws<std::domain_error>([&] { solve(make_small_mesh<T>(), T{ 2 }, {}, conditions, {}); }));
    };

    "nonlocal_stationary_solvers"_test = [solve] {
        const auto mesh = make_small_mesh<T>();
        const model_parameters<1, T> model{
            .influence = [](const T&, const T&) constexpr noexcept { return T{ 1 }; }, .local_weight = T{ 0.5 }
        };
        for (const bool is_variable : { false, true }) {
            const coefficient_t<T, 1> coefficient = is_variable
                ? coefficient_t<T, 1>{ spatial_dependency<T, 1>{
                      [](const point<T, 1>& x) constexpr noexcept { return T{ 1 } + x; } } }
                : coefficient_t<T, 1>{ T{ 2 } };
            const T left_load = is_variable ? T{ 0.5 } : T{ 0 };
            const T right_load = is_variable ? T{ 3.5 } : T{ 3 };
            const boundaries conditions{
                std::make_unique<Robin<T>>(T{ 3 }, left_load / T{ 3 }), std::make_unique<Load<T>>(right_load)
            };
            const solver_parameters parameters{
                .right_part = [is_variable](const T) constexpr noexcept { return is_variable ? T{ -0.5 } : T{ 0 }; }
            };
            // Use the K matrices from base/assemblers.cpp: c = 2 gives symmetric K; c = 1 + x gives nonsymmetric K.
            // For u = (1, 2, 3), (K + diag(3, 0, 0))*u is (0, 0, 3) or (1/4, -1/2, 13/4), respectively.
            static constexpr auto Expected = std::array{ T{ 1 }, T{ 2 }, T{ 3 } };
            check_values(solve(mesh, coefficient, model, conditions, parameters), Expected);
        }
    };
}

} // namespace nonlocal::unit_tests

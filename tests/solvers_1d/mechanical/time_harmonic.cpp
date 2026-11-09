#include <metamath/metamath.hpp>
#include <solvers/solver_1d/mechanical/harmonic_mechanical_equation_solver_1d.hpp>
#include <tests/solvers_1d/utils/solution_checks.hpp>

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"mechanical_time_harmonic_1d"> _ = [] {
    using namespace boost::ut;
    using namespace nonlocal;
    using namespace nonlocal::mesh;
    using namespace nonlocal::unit_tests;
    using namespace nonlocal::solver_1d::mechanical;
    using namespace metamath::constants;
    using namespace metamath::finite_element;
    using namespace metamath::functions;
    using T = double;
    using I = int64_t;

    static constexpr T Length = T{ 1 };
    static constexpr T Left_Stiffness = T{ 2 };
    static constexpr T Right_Stiffness = T{ 3 };
    static constexpr T Frequency = T{ 1e8 };
    static constexpr T young_modulus = Frequency * Frequency;
    static constexpr auto Expected_Displacement = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
    static constexpr auto Expected_Stress = [](const T x) { return young_modulus / ((6. - 3. * x) * std::cbrt(6. - 3. * x)); };
    // density * frequency^2 / youngs_modulus = 1
    const parameters_1d<T> parameters = { { .physical = { .youngs_modulus = young_modulus, .density = T{ 1 } } } };
    const time_harmonic_equation_parameters_1d<T> additional_parameters{
        .right_part =
            [&](const T x) {
                return -4. * young_modulus / (9. * power<2>(2. - x) * std::cbrt(6. - 3. * x)) -
                       Frequency * Frequency * Expected_Displacement(x);
            },
        .initial_distribution = [](const T x) { return Expected_Displacement(Length); },
        .frequency = Frequency
    };

    std::unordered_map<std::string, mechanical_boundaries_conditions_1d<T>> boundaries_conditions;
    boundaries_conditions["displacement_displacement"] = { std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
                                                           std::make_unique<displacement_1d<T>>(Expected_Displacement(Length)) };
    boundaries_conditions["displacement_force"] = { std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
                                                    std::make_unique<pressure_1d<T>>(Expected_Stress(Length)) };
    boundaries_conditions["displacement_spring"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };
    boundaries_conditions["force_displacement"] = { std::make_unique<pressure_1d<T>>(-Expected_Stress(0.)),
                                                    std::make_unique<displacement_1d<T>>(Expected_Displacement(Length)) };
    boundaries_conditions["force_spring"] = {
        std::make_unique<pressure_1d<T>>(-Expected_Stress(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };
    boundaries_conditions["spring_displacement"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
    boundaries_conditions["spring_force"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<pressure_1d<T>>(Expected_Stress(Length))
    };
    boundaries_conditions["spring_spring"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };

    for (const auto& [test_name, conditions] : boundaries_conditions) {
        boost::ut::test(test_name) = [&parameters, &conditions, &additional_parameters] {
            // Mass entries contain N_i*N_j, so retain quadrature order 2.
            static constexpr size_t Quadrature_Order = 2;
            check_mesh_convergence<T>(
                Length, Quadrature_Order,
                [&](const auto& mesh) {
                    return harmonic_mechanical_equation_solver_1d<T, I>(mesh, parameters, conditions, additional_parameters);
                },
                [](const auto& solution) constexpr noexcept -> const auto& { return solution.displacement(); },
                [](auto& solution) -> const auto& { return solution.calc_stress(); },
                Expected_Displacement, Expected_Stress);
        };
    }

    "symmetric_indefinite_operator"_test = [] {
        const auto mesh = std::make_shared<mesh_1d<T>>(make_element_1d_integrated<T>(1, 2),
            std::vector<segment_data<T>>{ { .length = T{ 1 }, .elements = 1 } });
        const parameters_1d<T> parameters{ { .physical = { .youngs_modulus = T{ 1 }, .density = T{ 3 } } } };
        // For omega = 1, K - M = [[0, -1.5], [-1.5, 0]]. It maps u = (1, 2) to f = (-3, -1.5).
        const mechanical_boundaries_conditions_1d<T> conditions{
            std::make_unique<pressure_1d<T>>(T{ -3 }), std::make_unique<pressure_1d<T>>(T{ -1.5 })
        };
        const auto solution = harmonic_mechanical_equation_solver_1d<T, I>(mesh, parameters, conditions, { .frequency = T{ 1 } });
        check_values(solution.displacement(), std::array{ T{ 1 }, T{ 2 } });
    };

    "nonlocal_operator_with_mass"_test = [] {
        const auto mesh = make_small_mesh<T>();
        const parameters_1d<T> parameters{ {
            .model = {
                .influence = [](const T&, const T&) constexpr noexcept { return T{ 1 }; },
                .local_weight = T{ 0.5 }
            },
            .physical = {
                .youngs_modulus = spatial_dependency<T, 1>{ [](const point<T, 1>& x) constexpr noexcept { return T{ 1 } + x; } },
                .density = T{ 6 }
            }
        } };
        // Use the nonsymmetric K from base/assemblers.cpp and the full symmetric M for rho = 6.
        // For omega = 1 and u = (1, 2, 3), (K - M)*u = (-6.75, -12.5, -4.75).
        const mechanical_boundaries_conditions_1d<T> conditions{
            std::make_unique<pressure_1d<T>>(T{ -0.5 }), std::make_unique<pressure_1d<T>>(T{ 1.5 })
        };
        const auto solution = harmonic_mechanical_equation_solver_1d<T, I>(
            mesh, parameters, conditions,
            { .right_part = [](const T) constexpr noexcept { return T{ -12.5 }; }, .frequency = T{ 1 } });
        check_values(solution.displacement(), std::array{ T{ 1 }, T{ 2 }, T{ 3 } });
    };
};

} // namespace

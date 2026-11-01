#include <metamath/metamath.hpp>
#include <solvers/solver_1d/thermal/stationary_heat_equation_solver_1d.hpp>
#include <tests/solvers_1d/utils/stationary_checks.hpp>

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"thermal_stationary_1d"> _ = [] {
    using namespace boost::ut;
    using namespace nonlocal;
    using namespace nonlocal::mesh;
    using namespace nonlocal::unit_tests;
    using namespace nonlocal::solver_1d::thermal;
    using namespace metamath::constants;
    using namespace metamath::finite_element;
    using namespace metamath::functions;
    using T = double;
    using I = int64_t;

    static constexpr T Length = T{ 1 };
    static constexpr T Left_Heat_Transfer = T{ 2 };
    static constexpr T Right_Heat_Transfer = T{ 3 };
    static constexpr auto Expected_Temperature = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
    static constexpr auto Expected_Flux = [](const T x) { return -1. / ((6. - 3. * x) * std::cbrt(6. - 3. * x)); };
    const parameters_1d<T> parameters = { { .physical = { .conductivity = T{ 1 } } } };
    const stationary_equation_parameters_1d<T> additional_parameters{
        .right_part = [](const T x) { return -4. / (9 * power<2>(2 - x) * std::cbrt(6. - 3. * x)); },
        .initial_distribution = [](const T x) { return Expected_Temperature(Length); }
    };

    std::unordered_map<std::string, thermal_boundaries_conditions_1d<T>> boundaries_conditions;
    boundaries_conditions["temperature_temperature"] = { std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
                                                         std::make_unique<temperature_1d<T>>(Expected_Temperature(Length)) };
    boundaries_conditions["temperature_flux"] = { std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
                                                  std::make_unique<flux_1d<T>>(-Expected_Flux(Length)) };
    boundaries_conditions["temperature_convection"] = { std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
                                                        std::make_unique<convection_1d<T>>(
                                                            Right_Heat_Transfer, -Expected_Flux(Length) / Right_Heat_Transfer +
                                                                                     Expected_Temperature(Length)) };
    boundaries_conditions["temperature_radiation"] = { std::make_unique<temperature_1d<T>>(Expected_Temperature(0)),
                                                       std::make_unique<radiation_1d<T>>(
                                                           -Expected_Flux(Length) / (Stefan_Boltzmann_Constant<T> *
                                                                                     power<4>(Expected_Temperature(Length)))) };
    boundaries_conditions["flux_temperature"] = { std::make_unique<flux_1d<T>>(Expected_Flux(0)),
                                                  std::make_unique<temperature_1d<T>>(Expected_Temperature(Length)) };
    boundaries_conditions["flux_convection"] = { std::make_unique<flux_1d<T>>(Expected_Flux(0)),
                                                 std::make_unique<convection_1d<T>>(Right_Heat_Transfer,
                                                                                    -Expected_Flux(Length) / Right_Heat_Transfer +
                                                                                        Expected_Temperature(Length)) };
    boundaries_conditions["flux_radiation"] = { std::make_unique<flux_1d<T>>(Expected_Flux(0)),
                                                std::make_unique<radiation_1d<T>>(
                                                    -Expected_Flux(Length) /
                                                    (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length)))) };
    boundaries_conditions["convection_temperature"] = {
        std::make_unique<convection_1d<T>>(Left_Heat_Transfer, Expected_Flux(0.) / Left_Heat_Transfer + Expected_Temperature(0.)),
        std::make_unique<temperature_1d<T>>(Expected_Temperature(Length))
    };
    boundaries_conditions["convection_flux"] = {
        std::make_unique<convection_1d<T>>(Left_Heat_Transfer, Expected_Flux(0.) / Left_Heat_Transfer + Expected_Temperature(0.)),
        std::make_unique<flux_1d<T>>(-Expected_Flux(Length))
    };
    boundaries_conditions["convection_convection"] = {
        std::make_unique<convection_1d<T>>(Left_Heat_Transfer, Expected_Flux(0.) / Left_Heat_Transfer + Expected_Temperature(0.)),
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer,
                                           -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["convection_radiation"] = {
        std::make_unique<convection_1d<T>>(Left_Heat_Transfer, Expected_Flux(0.) / Left_Heat_Transfer + Expected_Temperature(0.)),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) /
                                          (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
    };
    boundaries_conditions["radiation_temperature"] = {
        std::make_unique<radiation_1d<T>>(Expected_Flux(0) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(0)))),
        std::make_unique<temperature_1d<T>>(Expected_Temperature(Length))
    };
    boundaries_conditions["radiation_flux"] = {
        std::make_unique<radiation_1d<T>>(Expected_Flux(0) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(0)))),
        std::make_unique<flux_1d<T>>(-Expected_Flux(Length))
    };
    boundaries_conditions["radiation_convection"] = {
        std::make_unique<radiation_1d<T>>(Expected_Flux(0) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(0)))),
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer,
                                           -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["radiation_radiation"] = {
        std::make_unique<radiation_1d<T>>(Expected_Flux(0) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(0)))),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) /
                                          (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
    };

    for (const auto& [test_name, conditions] : boundaries_conditions) {
        boost::ut::test(test_name) = [&parameters, &conditions, &additional_parameters] {
            static constexpr size_t Quadrature_Order = 1;
            check_mesh_convergence<T>(
                Length, Quadrature_Order,
                [&](const auto& mesh) {
                    return stationary_heat_equation_solver_1d<T, I>(mesh, parameters, conditions, additional_parameters);
                },
                [](const auto& solution) constexpr noexcept -> const auto& { return solution.temperature(); },
                [](auto& solution) -> const auto& { return solution.calc_flux(); },
                Expected_Temperature, Expected_Flux);
        };
    }

    const auto solve_example = [](const std::shared_ptr<mesh_1d<T>>& mesh, const coefficient_t<T, 1>& coefficient,
                                  const model_parameters<1, T>& model, const thermal_boundaries_conditions_1d<T>& boundaries,
                                  const nonlocal::solver_1d::stationary_solver_parameters_1d<T>& settings) -> std::vector<T> {
        const parameters_1d<T> parameters{ { .model = model, .physical = { .conductivity = coefficient } } };
        const stationary_equation_parameters_1d<T> additional_parameters{
            .right_part = settings.right_part,
            .initial_distribution = settings.initial_distribution,
            .tolerance = settings.tolerance,
            .max_iterations = settings.max_iterations,
            .energy = settings.integral_value
        };
        return stationary_heat_equation_solver_1d<T, I>(mesh, parameters, boundaries, additional_parameters).temperature();
    };
    register_stationary_examples<T, physics_t::THERMAL, temperature_1d, flux_1d, convection_1d>(solve_example);
};

} // namespace

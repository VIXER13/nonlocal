#include "utils.hpp"

#include <metamath/metamath.hpp>
#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_1d/thermal/stationary_heat_equation_solver_1d.hpp>

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

    static constexpr T Length = T{1};
    static constexpr T Left_Heat_Transfer = T{2};
    static constexpr T Right_Heat_Transfer = T{3};
    static constexpr auto Expected_Temperature = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
    static constexpr auto Expected_Flux = [](const T x) { return -1. / ((6. - 3. * x) * std::cbrt(6. - 3. * x)); };
    const parameters_1d<T> parameters = {{ .physical = { .conductivity = T{1} } }};
    const stationary_equation_parameters_1d<T> additional_parameters{
        .right_part = [](const T x) {  return -4. / (9 * power<2>(2 - x) * std::cbrt(6. - 3. * x)); },
        .initial_distribution = [](const T x) { return Expected_Temperature(Length); }
    };
    
    std::unordered_map<std::string, thermal_boundaries_conditions_1d<T>> boundaries_conditions;
    boundaries_conditions["temperature_temperature"] = {
        std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
        std::make_unique<temperature_1d<T>>(Expected_Temperature(Length))
    };
    boundaries_conditions["temperature_flux"] = {
        std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
        std::make_unique<flux_1d<T>>(-Expected_Flux(Length))
    };
    boundaries_conditions["temperature_convection"] = {
        std::make_unique<temperature_1d<T>>(Expected_Temperature(0.)),
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer, -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["temperature_radiation"] = {
        std::make_unique<temperature_1d<T>>(Expected_Temperature(0)),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
    };
    boundaries_conditions["flux_temperature"] = {
        std::make_unique<flux_1d<T>>(Expected_Flux(0)),
        std::make_unique<temperature_1d<T>>(Expected_Temperature(Length))
    };
    boundaries_conditions["flux_convection"] = {
        std::make_unique<flux_1d<T>>(Expected_Flux(0)),
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer, -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["flux_radiation"] = {
        std::make_unique<flux_1d<T>>(Expected_Flux(0)),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
    };
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
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer, -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["convection_radiation"] = {
        std::make_unique<convection_1d<T>>(Left_Heat_Transfer, Expected_Flux(0.) / Left_Heat_Transfer + Expected_Temperature(0.)),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
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
        std::make_unique<convection_1d<T>>(Right_Heat_Transfer, -Expected_Flux(Length) / Right_Heat_Transfer + Expected_Temperature(Length))
    };
    boundaries_conditions["radiation_radiation"] = {
        std::make_unique<radiation_1d<T>>(Expected_Flux(0) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(0)))),
        std::make_unique<radiation_1d<T>>(-Expected_Flux(Length) / (Stefan_Boltzmann_Constant<T> * power<4>(Expected_Temperature(Length))))
    };

    for(const auto& node : boundaries_conditions) {
        const auto& test_name = node.first;
        const auto& conditions = node.second;
        boost::ut::test(test_name) = [&parameters, &conditions, &additional_parameters] {
            T prev_temperature_error = std::numeric_limits<T>::max();
            T prev_flux_error = std::numeric_limits<T>::max();
            for(const size_t elements : {10, 20, 40}) {
                static constexpr size_t Order = 1;
                using quadrature = quadrature_1d<T, gauss, Order>;
                using element_1d = element_1d_integrate<T, lagrangian_element_1d, Order>;
                const auto mesh = std::make_shared<mesh_1d<T>>(
                    std::make_unique<element_1d>(quadrature{}),
                    std::vector<segment_data<T>>{{ .length = Length, .elements = elements }});
                auto solution = stationary_heat_equation_solver_1d<T, I>(mesh, parameters, conditions, additional_parameters);
                solution.calc_flux();
                const auto expected_temperature_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Temperature);
                const auto expected_flux_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Flux);
                const T temperature_error = max_error(solution.temperature(), expected_temperature_discrete) / max_norm(expected_temperature_discrete);
                const T flux_error = max_error(solution.flux(), expected_flux_discrete) / max_norm(expected_flux_discrete);
                expect(lt(temperature_error, prev_temperature_error));
                expect(lt(flux_error, prev_flux_error));
                prev_temperature_error = temperature_error;
                prev_flux_error = flux_error;
            }
        };
    }
};
    
}
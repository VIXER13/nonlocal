#include <metamath/metamath.hpp>
#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_1d/mechanical/harmonic_mechanical_equation_solver_1d.hpp>
#include <tests/utils/error.hpp>

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

    static constexpr T Length = T{1};
    static constexpr T Left_Stiffness = T{2};
    static constexpr T Right_Stiffness = T{3};
    static constexpr T Frequency = T{1e8};
    static constexpr T young_modulus = Frequency * Frequency;
    static constexpr auto Expected_Displacement = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
    static constexpr auto Expected_Stress = [](const T x) { return young_modulus / ((6. - 3. * x) * std::cbrt(6. - 3. * x)); };
    // density * frequency^2 / youngs_modulus = 1
    const parameters_1d<T> parameters = {{ .physical = { .youngs_modulus = young_modulus, .density = T{1} } }};
    const time_harmonic_equation_parameters_1d<T> additional_parameters{
        .right_part = [&](const T x) {  return -4. * young_modulus / (9. * power<2>(2. - x) * std::cbrt(6. - 3. * x)) - Frequency * Frequency * Expected_Displacement(x); },
        .initial_distribution = [](const T x) { return Expected_Displacement(Length); },
        .frequency = Frequency
    };
    
    std::unordered_map<std::string, mechanical_boundaries_conditions_1d<T>> boundaries_conditions;
    boundaries_conditions["displacement_displacement"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
    boundaries_conditions["displacement_force"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<pressure_1d<T>>(Expected_Stress(Length))
    };
    boundaries_conditions["displacement_spring"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };
    boundaries_conditions["force_displacement"] = {
        std::make_unique<pressure_1d<T>>(-Expected_Stress(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
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

    for(const auto& [test_name, conditions] : boundaries_conditions) {
        boost::ut::test(test_name) = [&parameters, &conditions, &additional_parameters] {
            constexpr size_t number_of_metrics = 4;
            std::array<T, number_of_metrics> prev_errors{};
            std::array<T, number_of_metrics> curr_errors{};
            std::fill(prev_errors.begin(), prev_errors.end(), std::numeric_limits<T>::max());
            for(const size_t elements : {10, 20, 40}) {
                static constexpr size_t Element_Order = 1;
                // Mass matrix must be intergrated with quadratures one order higher then Stiffness
                // As mesh contains one set of quadrature points for both matrices one need to increase it's order for Mass matrix
                static constexpr size_t Quadrature_Order = 2;
                const auto mesh = std::make_shared<mesh_1d<T>>(
                    metamath::finite_element::make_element_1d_integrated<T>(Element_Order, Quadrature_Order),
                    std::vector<segment_data<T>>{{ .length = Length, .elements = elements }});
                auto solution = harmonic_mechanical_equation_solver_1d<T, I>(mesh, parameters, conditions, additional_parameters);
                solution.calc_stress();
                const auto Expected_Displacement_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Displacement);
                const auto Expected_Stress_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Stress);
                curr_errors[0] = L2_norm  (solution.displacement(), Expected_Displacement_discrete);
                curr_errors[1] = L2_norm  (solution.stress(),       Expected_Stress_discrete);
                curr_errors[2] = max_error(solution.displacement(), Expected_Displacement_discrete);
                curr_errors[3] = max_error(solution.stress(),       Expected_Stress_discrete);
                for (size_t i = 0; i < number_of_metrics; ++i)
                    expect(lt(curr_errors[i], prev_errors[i]));
                std::swap(curr_errors, prev_errors);
            }
        };
    }
};
    
}
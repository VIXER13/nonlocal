
#include <metamath/metamath.hpp>
#include <mesh/mesh_1d/mesh_1d_utils.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_1d/mechanical/stationary_mechanical_equation_solver_1d.hpp>
#include <tests/utils/error.hpp>

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"mechanical_stationary_1d"> _ = [] {
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
    static constexpr auto Expected_Displacement = [](const T x) { return 1. / std::cbrt(6. - 3. * x); };
    static constexpr auto Expected_Stress = [](const T x) { return 1. / ((6. - 3. * x) * std::cbrt(6. - 3. * x)); };
    const parameters_1d<T> parameters = {{ .physical = { .youngs_modulus = T{1} } }};
    const stationary_equation_parameters_1d<T> additional_parameters{
        .right_part = [](const T x) {  return -4. / (9 * power<2>(2 - x) * std::cbrt(6. - 3. * x)); },
        .initial_distribution = [](const T x) { return Expected_Displacement(Length); }
    };
    
    std::unordered_map<std::string, mechanical_boundaries_conditions_1d<T>> boundaries_conditions;
    boundaries_conditions["displacement_displacement"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
    boundaries_conditions["displacement_force"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<normal_force_1d<T>>(Expected_Stress(Length))
    };
    boundaries_conditions["displacement_spring"] = {
        std::make_unique<displacement_1d<T>>(Expected_Displacement(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };
    boundaries_conditions["force_displacement"] = {
        std::make_unique<normal_force_1d<T>>(-Expected_Stress(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
    boundaries_conditions["force_spring"] = {
        std::make_unique<normal_force_1d<T>>(-Expected_Stress(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };
    boundaries_conditions["spring_displacement"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<displacement_1d<T>>(Expected_Displacement(Length))
    };
    boundaries_conditions["spring_force"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<normal_force_1d<T>>(Expected_Stress(Length))
    };
    boundaries_conditions["spring_spring"] = {
        std::make_unique<spring_1d<T>>(Left_Stiffness, -Expected_Stress(0.) / Left_Stiffness + Expected_Displacement(0.)),
        std::make_unique<spring_1d<T>>(Right_Stiffness, Expected_Stress(Length) / Right_Stiffness + Expected_Displacement(Length))
    };

    for(const auto& node : boundaries_conditions) {
        const auto& test_name = node.first;
        const auto& conditions = node.second;
        boost::ut::test(test_name) = [&parameters, &conditions, &additional_parameters] {
            T prev_displacement_error = std::numeric_limits<T>::max();
            T prev_stress_error = std::numeric_limits<T>::max();
            for(const size_t elements : {10, 20, 40}) {
                static constexpr size_t Order = 1;
                using quadrature = quadrature_1d<T, gauss, Order>;
                using element_integrate_1d = element_1d_integrate<T>;
                const auto mesh = std::make_shared<mesh_1d<T>>(
                    std::make_unique<element_integrate_1d>(
                        std::make_unique<element_1d<T, lagrangian_element_1d, Order>>(),
                        quadrature{}),
                    std::vector<segment_data<T>>{{ .length = Length, .elements = elements }});
                auto solution = stationary_mechanical_equation_solver_1d<T, I>(mesh, parameters, conditions, additional_parameters);
                solution.calc_stress();
                const auto Expected_Displacement_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Displacement);
                const auto Expected_Stress_discrete = nonlocal::mesh::utils::discrete<T>(*mesh, Expected_Stress);
                const T displacement_error = L2_norm(solution.displacement(), Expected_Displacement_discrete);
                const T stress_error = L2_norm(solution.stress(), Expected_Stress_discrete);
                expect(lt(displacement_error, prev_displacement_error));
                expect(lt(stress_error, prev_stress_error));
                prev_displacement_error = displacement_error;
                prev_stress_error = stress_error;
            }
        };
    }
};
    
}
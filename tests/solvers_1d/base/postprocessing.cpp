#include <solvers/solver_1d/mechanical/mechanical_equation_solution_1d.hpp>
#include <solvers/solver_1d/thermal/heat_equation_solution_1d.hpp>
#include <tests/solvers_1d/utils/solution_checks.hpp>

#include <boost/ut.hpp>

#include <cmath>

namespace {

using T = double;
using namespace boost::ut;
using namespace nonlocal;
using namespace nonlocal::solver_1d;
using namespace nonlocal::unit_tests;

const suite<"postprocessing_1d"> _ = [] {
    "weighted_gradient_with_two_segments"_test = [] {
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(
            metamath::finite_element::make_element_1d_integrated<T>(1, 1),
            std::vector<mesh::segment_data<T>>{ { .length = T{ 1 }, .elements = 1 },
                                               { .length = T{ 4 }, .search_radius = T{ 2 }, .elements = 2 } });
        const std::vector<T> values{ T{ 0 }, T{ 1 }, T{ 3 }, T{ 7 } };
        const model_parameters<1, T> nonlocal_model{
            .influence = [](const T, const T) constexpr noexcept { return T{ 0.25 }; }, .local_weight = T{ 0.25 }
        };
        const solution_dependency<T, 1> coefficient = [](const point<T, 1>&, const T value) constexpr noexcept {
            return value;
        };
        const mechanical::parameters_1d<T> mechanical_parameters{
            { .physical = { .youngs_modulus = T{ 2 }, .density = T{ 1 } } },
            { .model = nonlocal_model, .physical = { .youngs_modulus = coefficient, .density = T{ 1 } } }
        };
        const thermal::parameters_1d<T> thermal_parameters{
            { .physical = { .conductivity = T{ 2 } } },
            { .model = nonlocal_model, .physical = { .conductivity = coefficient } }
        };
        mechanical::mechanical_equation_solution_1d<T> mechanical_solution{ mesh, mechanical_parameters, values };
        thermal::heat_equation_solution_1d<T> thermal_solution{ mesh, thermal_parameters, values };

        // The local weighted gradients are 2, 2 and 10. On the second segment,
        // the nonlocal average is 6, so the mixed values are 5 and 7.
        // At shared nodes we average the adjacent element values.
        static constexpr std::array Expected_Stress{ T{ 2 }, T{ 3.5 }, T{ 6 }, T{ 7 } };
        static constexpr std::array Expected_Flux{ T{ -2 }, T{ -3.5 }, T{ -6 }, T{ -7 } };
        check_values(mechanical_solution.calc_stress(), Expected_Stress);
        check_values(thermal_solution.calc_flux(), Expected_Flux);
    };

    "heat_flux_uses_previous_step"_test = [] {
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(
            metamath::finite_element::make_element_1d_integrated<T>(1, 1),
            std::vector<mesh::segment_data<T>>{ { .length = T{ 1 }, .elements = 1 } });
        const thermal::parameters_1d<T> parameters{
            { .physical = { .conductivity = T{ 2 }, .relaxation_time = T{ 1 } } }
        };
        const T time_step = std::log(T{ 2 });
        thermal::heat_equation_solution_1d<T> solution{ mesh, parameters, std::vector{ T{ 0 }, T{ 1 } }, time_step };

        // The relaxation factor halves the conductivity. The initial flux is -1.
        check_values(solution.calc_flux(), std::array{ T{ -1 }, T{ -1 } });
        solution.temperature(std::vector{ T{ 0 }, T{ 2 } });
        solution.time(time_step);
        // The next flux is -2 plus half the previous flux.
        check_values(solution.calc_flux(), std::array{ T{ -2.5 }, T{ -2.5 } });
    };
};

} // namespace

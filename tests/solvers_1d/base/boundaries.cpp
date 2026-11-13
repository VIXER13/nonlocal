#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_nonlinear_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_second_kind_1d.hpp>
#include <solvers/solver_1d/mechanical/init_problem_settings.hpp>
#include <solvers/solver_1d/mechanical/spring_condition_1d.hpp>
#include <solvers/solver_1d/thermal/convection_condition_1d.hpp>
#include <solvers/solver_1d/thermal/init_problem_settings.hpp>
#include <solvers/solver_1d/thermal/radiation_condition_1d.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using I = int64_t;
using vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using sparse_matrix = Eigen::SparseMatrix<T, Eigen::RowMajor, I>;
using namespace nonlocal;
using namespace nonlocal::solver_1d;
using namespace boost::ut;

static constexpr T Epsilon = T{ 1e-14 };
static constexpr T Sigma = metamath::constants::Stefan_Boltzmann_Constant<T>;

sparse_matrix make_matrix() {
    sparse_matrix matrix{ 3, 3 };
    static constexpr auto Diagonal = std::array{ T{ 10 }, T{ 20 }, T{ 30 } };
    for (const size_t i : std::ranges::iota_view{ 0zu, Diagonal.size() })
        matrix.insert(i, i) = Diagonal[i];
    matrix.makeCompressed();
    return matrix;
}

template<class Vector, size_t N>
void check_values(const Vector& values, const std::array<T, N>& expected) {
    for (const size_t i : std::ranges::iota_view{ 0zu, N })
        expect(approx(values[i], expected[i], Epsilon));
}

void check_diagonal(const sparse_matrix& matrix, const std::array<T, 3>& expected) {
    for (const size_t i : std::ranges::iota_view{ 0zu, 3zu })
        for (const size_t j : std::ranges::iota_view{ 0zu, 3zu })
            expect(approx(matrix.coeff(i, j), i == j ? expected[i] : T{ 0 }, Epsilon));
}

const suite<"boundaries_1d"> _ = [] {
    "prescribed_values_and_rhs"_test = [] {
        const thermal::thermal_boundaries_conditions_1d<T> thermal_boundaries{
            std::make_unique<thermal::temperature_1d<T>>(T{ 3 }),
            std::make_unique<thermal::temperature_1d<T>>(T{ 5 })
        };
        const mechanical::mechanical_boundaries_conditions_1d<T> mechanical_boundaries{
            std::make_unique<mechanical::displacement_1d<T>>(T{ 3 }),
            std::make_unique<mechanical::displacement_1d<T>>(T{ 5 })
        };
        const std::array<std::unordered_map<size_t, T>, 2> bound{
            std::unordered_map<size_t, T>{ { 1zu, T{ -2 } } }, std::unordered_map<size_t, T>{ { 1zu, T{ -4 } } }
        };
        std::array<T, 3> thermal_rhs{ 0., 1., 0. };
        auto mechanical_rhs = thermal_rhs;
        boundary_condition_first_kind_1d(thermal_rhs, bound, thermal_boundaries);
        boundary_condition_first_kind_1d(mechanical_rhs, bound, mechanical_boundaries);
        // u(0) = 3, u(2) = 5: the middle load becomes 1 - (-2)*3 - (-4)*5 = 27.
        static constexpr auto Expected = std::array{ T{ 3 }, T{ 27 }, T{ 5 } };
        check_values(thermal_rhs, Expected);
        check_values(mechanical_rhs, Expected);
    };

    "loads_and_neumann_index"_test = [] {
        const thermal::thermal_boundaries_conditions_1d<T> thermal_boundaries{
            std::make_unique<thermal::flux_1d<T>>(T{ 2 }), std::make_unique<thermal::flux_1d<T>>(T{ -3 })
        };
        const mechanical::mechanical_boundaries_conditions_1d<T> mechanical_boundaries{
            std::make_unique<mechanical::pressure_1d<T>>(T{ 2 }), std::make_unique<mechanical::pressure_1d<T>>(T{ -3 })
        };
        std::array<T, 4> thermal_rhs{};
        std::array<T, 4> mechanical_rhs{};
        boundary_condition_second_kind_1d(thermal_rhs, thermal_boundaries, true);
        boundary_condition_second_kind_1d(mechanical_rhs, mechanical_boundaries, true);
        // Three physical nodes and one multiplier: the right load belongs to node 2, not the multiplier.
        static constexpr auto Expected = std::array{ T{ 2 }, T{ 0 }, T{ -3 }, T{ 0 } };
        check_values(thermal_rhs, Expected);
        check_values(mechanical_rhs, Expected);
        const thermal::parameters_1d<T> thermal_parameters;
        const mechanical::parameters_1d<T> mechanical_parameters;
        expect(thermal::init_problem_settings(thermal_parameters, thermal_boundaries, true).is_neumann);
        expect(!thermal::init_problem_settings(thermal_parameters, thermal_boundaries, false).is_neumann);
        expect(mechanical::init_problem_settings(mechanical_parameters, mechanical_boundaries, true).is_neumann);
        expect(!mechanical::init_problem_settings(mechanical_parameters, mechanical_boundaries, false).is_neumann);
    };

    "convection_and_spring"_test = [] {
        const thermal::thermal_boundaries_conditions_1d<T> thermal_boundaries{
            std::make_unique<thermal::convection_1d<T>>(T{ 3 }, T{ 10 }),
            std::make_unique<thermal::flux_1d<T>>(T{ 5 })
        };
        const mechanical::mechanical_boundaries_conditions_1d<T> mechanical_boundaries{
            std::make_unique<mechanical::spring_1d<T>>(T{ 3 }, T{ 10 }),
            std::make_unique<mechanical::pressure_1d<T>>(T{ 5 })
        };
        auto thermal_matrix = make_matrix();
        auto mechanical_matrix = make_matrix();
        std::array<T, 3> thermal_rhs{};
        std::array<T, 3> mechanical_rhs{};
        thermal::convection_condition_1d(thermal_matrix, thermal_boundaries);
        mechanical::spring_condition_1d(mechanical_matrix, mechanical_boundaries);
        boundary_condition_second_kind_1d(thermal_rhs, thermal_boundaries);
        boundary_condition_second_kind_1d(mechanical_rhs, mechanical_boundaries);
        // Robin coefficient 3 adds 3 to A(0,0) and 3*10 to f(0); the right load is 5.
        static constexpr auto Expected_Diagonal = std::array{ T{ 13 }, T{ 20 }, T{ 30 } };
        static constexpr auto Expected_Load = std::array{ T{ 30 }, T{ 0 }, T{ 5 } };
        check_diagonal(thermal_matrix, Expected_Diagonal);
        check_diagonal(mechanical_matrix, Expected_Diagonal);
        check_values(thermal_rhs, Expected_Load);
        check_values(mechanical_rhs, Expected_Load);
    };

    "combined_conditions"_test = [] {
        const thermal::combined_flux_1d<T> thermal_condition{ T{ 2 }, T{ 3 }, T{ 10 }, T{ 0.5 } };
        const mechanical::combined_loading_1d<T> mechanical_condition{ T{ 2 }, T{ 3 }, T{ 10 } };
        // Both linear loads are 2 + 3*10 = 32; radiation is a separate nonlinear term.
        const thermal::thermal_boundary_condition_1d<T>& thermal_boundary = thermal_condition;
        const mechanical::mechanical_boundary_condition_1d<T>& mechanical_boundary = mechanical_condition;
        expect(eq(thermal_boundary(), T{ 32 }));
        expect(eq(mechanical_boundary(), T{ 32 }));
        const auto* const radiation = dynamic_cast<const thermal::radiation_1d<T>*>(&thermal_boundary);
        expect(radiation != nullptr) << fatal;
        expect(eq(radiation->emissivity(), T{ 0.5 }));
        std::array<T, 3> thermal_rhs{};
        auto mechanical_rhs = thermal_rhs;
        auto thermal_matrix = make_matrix();
        auto mechanical_matrix = make_matrix();
        boundary_condition_second_kind_1d(thermal_rhs, thermal_boundary, 0);
        boundary_condition_second_kind_1d(mechanical_rhs, mechanical_boundary, 0);
        thermal::convection_condition_1d(thermal_matrix, thermal_boundary, 0);
        mechanical::spring_condition_1d(mechanical_matrix, mechanical_boundary, 0);
        static constexpr auto Expected_Load = std::array{ T{ 32 }, T{ 0 }, T{ 0 } };
        static constexpr auto Expected_Diagonal = std::array{ T{ 13 }, T{ 20 }, T{ 30 } };
        check_values(thermal_rhs, Expected_Load);
        check_values(mechanical_rhs, Expected_Load);
        check_diagonal(thermal_matrix, Expected_Diagonal);
        check_diagonal(mechanical_matrix, Expected_Diagonal);
    };

    "zero_values_keep_boundary_kinds"_test = [] {
        const thermal::parameters_1d<T> thermal_parameters;
        const mechanical::parameters_1d<T> mechanical_parameters;
        const thermal::thermal_boundaries_conditions_1d<T> prescribed{
            std::make_unique<thermal::temperature_1d<T>>(T{ 0 }), std::make_unique<thermal::flux_1d<T>>(T{ 0 })
        };
        const auto prescribed_settings = thermal::init_problem_settings(thermal_parameters, prescribed, true);
        expect(prescribed_settings.is_first_kind.front());
        expect(!prescribed_settings.is_neumann);
        // Zero coefficients retain the original condition kinds and solver branch selection.
        const thermal::thermal_boundaries_conditions_1d<T> thermal_boundaries{
            std::make_unique<thermal::combined_flux_1d<T>>(T{ 2 }, T{ 0 }, T{ 0 }, T{ 0 }),
            std::make_unique<thermal::flux_1d<T>>(T{ 0 })
        };
        const mechanical::mechanical_boundaries_conditions_1d<T> mechanical_boundaries{
            std::make_unique<mechanical::combined_loading_1d<T>>(T{ 2 }, T{ 0 }, T{ 0 }),
            std::make_unique<mechanical::pressure_1d<T>>(T{ 0 })
        };
        const auto thermal_settings = thermal::init_problem_settings(thermal_parameters, thermal_boundaries, true);
        expect(thermal_settings.is_nonlinear_boundary);
        expect(!thermal_settings.is_neumann);
        expect(!mechanical::init_problem_settings(mechanical_parameters, mechanical_boundaries, true).is_neumann);
    };

    "nonlinear_stationary_boundary"_test = [] {
        auto matrix = make_matrix();
        vector residual(3);
        residual << 5., 6., 7.;
        // g(u) = u^2 at u = 2: A(0,0) = 10 - g'(2) = 6, residual(0) = 5 - g(2) = 1.
        boundary_condition_nonlinear_1d<true>(matrix, residual, 0, T{ 2 }, T{ 4 }, T{ 4 });
        static constexpr auto Expected_Diagonal = std::array{ T{ 6 }, T{ 20 }, T{ 30 } };
        static constexpr auto Expected_Residual = std::array{ T{ 1 }, T{ 6 }, T{ 7 } };
        check_diagonal(matrix, Expected_Diagonal);
        check_values(residual, Expected_Residual);
    };

    "nonlinear_transient_boundary"_test = [] {
        auto matrix = make_matrix();
        vector rhs(3);
        rhs << 5., 6., 7.;
        // dt = 1/4: A(0,0) = 10 + g'(2)/4 = 11, f(0) = 5 + (g'(2)*2 - g(2))/4 = 6.
        boundary_condition_nonlinear_1d<false>(matrix, rhs, 0, T{ 2 }, T{ 4 }, T{ 4 }, T{ 0.25 });
        static constexpr auto Expected_Diagonal = std::array{ T{ 11 }, T{ 20 }, T{ 30 } };
        static constexpr auto Expected_Load = std::array{ T{ 6 }, T{ 6 }, T{ 7 } };
        check_diagonal(matrix, Expected_Diagonal);
        check_values(rhs, Expected_Load);
    };

    "radiation_in_both_modes"_test = [] {
        const thermal::thermal_boundaries_conditions_1d<T> boundaries{
            std::make_unique<thermal::radiation_1d<T>>(T{ 0.5 }), std::make_unique<thermal::flux_1d<T>>(T{ 0 })
        };
        // g(u) = 0.5*sigma*u^4: g(2) = 8*sigma and g'(2) = 16*sigma.
        auto stationary_matrix = make_matrix();
        auto transient_matrix = make_matrix();
        vector residual = vector::Constant(3, T{ 5 });
        vector rhs = residual;
        const vector previous = vector::Constant(3, T{ 2 });
        thermal::radiation_condition_1d<true>(stationary_matrix, residual, boundaries, previous);
        thermal::radiation_condition_1d<false>(transient_matrix, rhs, boundaries, previous, T{ 0.25 });
        static constexpr auto Stationary_Diagonal = std::array{ T{ 10 } - T{ 16 } * Sigma, T{ 20 }, T{ 30 } };
        static constexpr auto Transient_Diagonal = std::array{ T{ 10 } + T{ 4 } * Sigma, T{ 20 }, T{ 30 } };
        static constexpr auto Stationary_Residual = std::array{ T{ 5 } - T{ 8 } * Sigma, T{ 5 }, T{ 5 } };
        static constexpr auto Transient_Load = std::array{ T{ 5 } + T{ 6 } * Sigma, T{ 5 }, T{ 5 } };
        check_diagonal(stationary_matrix, Stationary_Diagonal);
        check_diagonal(transient_matrix, Transient_Diagonal);
        check_values(residual, Stationary_Residual);
        check_values(rhs, Transient_Load);
        const thermal::parameters_1d<T> parameters{ { .physical = { .conductivity = T{ 1 } } } };
        const auto settings = thermal::init_problem_settings(parameters, boundaries, true);
        expect(settings.is_nonlinear_boundary);
        expect(!settings.is_neumann);
    };
};

} // namespace

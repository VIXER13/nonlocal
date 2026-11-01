#include <solvers/solver_1d/base/boundary_condition_first_kind_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_nonlinear_1d.hpp>
#include <solvers/solver_1d/base/boundary_condition_robin_1d.hpp>
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

boundaries_data_1d<T> quadratic_boundary() {
    return { boundary_data_1d<T>{
                 .nonlinear = nonlinear_boundary_operator_1d<T>{
                     .value = [](const T u) constexpr noexcept { return u * u; },
                     .derivative = [](const T u) constexpr noexcept { return T{ 2 } * u; }
                 }
             },
             boundary_data_1d<T>{} };
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
        boundary_condition_first_kind_1d(thermal_rhs, bound, normalize_boundaries_1d(thermal_boundaries));
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
        const auto data = normalize_boundaries_1d(thermal_boundaries);
        std::array<T, 4> thermal_rhs{};
        std::array<T, 4> mechanical_rhs{};
        boundary_condition_second_kind_1d(thermal_rhs, data, true);
        boundary_condition_second_kind_1d(mechanical_rhs, mechanical_boundaries, true);
        // Three physical nodes and one multiplier: the right load belongs to node 2, not the multiplier.
        static constexpr auto Expected = std::array{ T{ 2 }, T{ 0 }, T{ -3 }, T{ 0 } };
        check_values(thermal_rhs, Expected);
        check_values(mechanical_rhs, Expected);
        expect(init_boundary_settings(problem_settings{}, data, true).is_neumann);
        expect(!init_boundary_settings(problem_settings{}, data, false).is_neumann);
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
        const auto thermal_data = thermal_condition.data();
        const auto mechanical_data = mechanical_condition.data();
        // Both linear loads are 2 + 3*10 = 32; radiation is a separate nonlinear term.
        expect(eq(thermal_data.load, T{ 32 }));
        expect(eq(mechanical_data.load, T{ 32 }));
        expect(thermal_data.robin_coefficient.has_value()) << fatal;
        expect(mechanical_data.robin_coefficient.has_value()) << fatal;
        expect(eq(*thermal_data.robin_coefficient, T{ 3 }));
        expect(eq(*mechanical_data.robin_coefficient, T{ 3 }));
        expect(thermal_data.nonlinear.has_value()) << fatal;
        expect(!mechanical_data.nonlinear);
        expect(approx(thermal_data.nonlinear->value(T{ 2 }), T{ 8 } * Sigma, Epsilon));
        std::array<T, 3> rhs{};
        auto matrix = make_matrix();
        boundary_condition_second_kind_1d(rhs, thermal_data, 0);
        boundary_condition_robin_1d(matrix, thermal_data, 0);
        static constexpr auto Expected_Load = std::array{ T{ 32 }, T{ 0 }, T{ 0 } };
        static constexpr auto Expected_Diagonal = std::array{ T{ 13 }, T{ 20 }, T{ 30 } };
        check_values(rhs, Expected_Load);
        check_diagonal(matrix, Expected_Diagonal);
    };

    "zero_values_keep_boundary_kinds"_test = [] {
        const boundaries_data_1d<T> prescribed{ boundary_data_1d<T>{ .prescribed_value = T{ 0 } }, boundary_data_1d<T>{} };
        const auto prescribed_settings = init_boundary_settings(problem_settings{}, prescribed, true);
        expect(prescribed_settings.is_first_kind.front());
        expect(!prescribed_settings.is_neumann);
        const thermal::combined_flux_1d<T> thermal_condition{ T{ 2 }, T{ 0 }, T{ 0 }, T{ 0 } };
        const mechanical::combined_loading_1d<T> mechanical_condition{ T{ 2 }, T{ 0 }, T{ 0 } };
        // Zero coefficients retain the original condition kinds and solver branch selection.
        const auto thermal_data = thermal_condition.data();
        const auto mechanical_data = mechanical_condition.data();
        expect(thermal_data.robin_coefficient.has_value());
        expect(thermal_data.nonlinear.has_value());
        expect(!thermal_data.is_pure_load());
        expect(mechanical_data.robin_coefficient.has_value());
        expect(!mechanical_data.is_pure_load());
    };

    "nonlinear_stationary_boundary"_test = [] {
        auto matrix = make_matrix();
        vector residual(3);
        residual << 5., 6., 7.;
        const vector previous = vector::Constant(3, T{ 2 });
        boundary_condition_nonlinear_1d<true>(matrix, residual, quadratic_boundary(), previous);
        // g(u) = u^2 at u = 2: A(0,0) = 10 - g'(2) = 6, residual(0) = 5 - g(2) = 1.
        static constexpr auto Expected_Diagonal = std::array{ T{ 6 }, T{ 20 }, T{ 30 } };
        static constexpr auto Expected_Residual = std::array{ T{ 1 }, T{ 6 }, T{ 7 } };
        check_diagonal(matrix, Expected_Diagonal);
        check_values(residual, Expected_Residual);
    };

    "nonlinear_transient_boundary"_test = [] {
        auto matrix = make_matrix();
        vector rhs(3);
        rhs << 5., 6., 7.;
        const vector previous = vector::Constant(3, T{ 2 });
        boundary_condition_nonlinear_1d<false>(matrix, rhs, quadratic_boundary(), previous, T{ 0.25 });
        // dt = 1/4: A(0,0) = 10 + g'(2)/4 = 11, f(0) = 5 + (g'(2)*2 - g(2))/4 = 6.
        static constexpr auto Expected_Diagonal = std::array{ T{ 11 }, T{ 20 }, T{ 30 } };
        static constexpr auto Expected_Load = std::array{ T{ 6 }, T{ 6 }, T{ 7 } };
        check_diagonal(matrix, Expected_Diagonal);
        check_values(rhs, Expected_Load);
    };

    "radiation_normalization_and_modes"_test = [] {
        const thermal::thermal_boundaries_conditions_1d<T> boundaries{
            std::make_unique<thermal::radiation_1d<T>>(T{ 0.5 }), std::make_unique<thermal::flux_1d<T>>(T{ 0 })
        };
        const auto data = normalize_boundaries_1d(boundaries);
        expect(data.front().nonlinear.has_value()) << fatal;
        // g(u) = 0.5*sigma*u^4: g(2) = 8*sigma and g'(2) = 16*sigma.
        expect(approx(data.front().nonlinear->value(T{ 2 }), T{ 8 } * Sigma, Epsilon));
        expect(approx(data.front().nonlinear->derivative(T{ 2 }), T{ 16 } * Sigma, Epsilon));
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
        const auto settings = thermal::init_problem_settings(parameters, data, true);
        expect(settings.is_nonlinear_boundary);
        expect(!settings.is_neumann);
    };
};

} // namespace

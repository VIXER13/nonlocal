#include <solvers/solver_1d/base/assemble_matrix_portrait.hpp>
#include <solvers/solver_1d/mechanical/mass_matrix_assembler.hpp>
#include <solvers/solver_1d/mechanical/stiffness_matrix_assembler.hpp>
#include <solvers/solver_1d/thermal/heat_capacity_assembler.hpp>
#include <solvers/solver_1d/thermal/thermal_conductivity_assembler.hpp>

#include <boost/ut.hpp>

namespace {

using T = double;
using I = int64_t;
using namespace nonlocal;
using namespace nonlocal::solver_1d;
using namespace boost::ut;

using matrix_3 = std::array<std::array<T, 3>, 3>;

std::shared_ptr<mesh::mesh_1d<T>> make_mesh() {
    return std::make_shared<mesh::mesh_1d<T>>(
        metamath::finite_element::make_element_1d_integrated<T>(1, 2),
        std::vector<mesh::segment_data<T>>{ { .length = T{ 2 }, .search_radius = T{ 2 }, .elements = 2 } });
}

// Expected arrays describe inner storage; for symmetric matrices, zeros below the diagonal mean omitted entries.
template<size_t N>
void check_matrix(const finite_element_matrix_1d<T, I>& matrix, const std::array<std::array<T, N>, N>& expected) {
    static constexpr T Epsilon = T{ 1e-14 };
    expect(eq(matrix.inner.rows(), I{ N })) << fatal;
    expect(eq(matrix.inner.cols(), I{ N })) << fatal;
    for (const size_t i : std::ranges::iota_view{ 0zu, N })
        for (const size_t j : std::ranges::iota_view{ 0zu, N })
            expect(approx(matrix.inner.coeff(i, j), expected[i][j], Epsilon)) << " at " << i << ", " << j;
}

template<class Assembler, class Parameters>
finite_element_matrix_1d<T, I> assemble_spatial(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                               const Parameters& parameters, const problem_settings& settings,
                                               const std::optional<std::vector<T>>& solution = std::nullopt) {
    finite_element_matrix_1d<T, I> matrix;
    init_matrix_portrait(matrix.inner, *mesh, settings);
    Assembler{ matrix, mesh }.calc_matrix(parameters, settings, solution);
    return matrix;
}

template<class Assembler, class Parameters>
finite_element_matrix_1d<T, I> assemble_mass(const std::shared_ptr<mesh::mesh_1d<T>>& mesh,
                                            const Parameters& parameters, const std::array<bool, 2>& is_first_kind) {
    finite_element_matrix_1d<T, I> matrix;
    const problem_settings settings{ .theories = { theory_t::LOCAL }, .is_first_kind = is_first_kind };
    init_matrix_portrait(matrix.inner, *mesh, settings);
    Assembler{ matrix, mesh }.calc_matrix(parameters, is_first_kind);
    return matrix;
}

void check_spatial(const std::shared_ptr<mesh::mesh_1d<T>>& mesh, const coefficient_t<T, 1>& coefficient,
                    const model_parameters<1, T>& model, const problem_settings& settings, const matrix_3& expected,
                    const std::optional<std::vector<T>>& solution = std::nullopt) {
    const thermal::parameters_1d<T> thermal_parameters{ { .model = model, .physical = { .conductivity = coefficient } } };
    const mechanical::parameters_1d<T> mechanical_parameters{
        { .model = model, .physical = { .youngs_modulus = coefficient, .density = T{ 1 } } }
    };
    check_matrix(assemble_spatial<thermal::thermal_conductivity_assembler_1d<T, I>>(
                     mesh, thermal_parameters, settings, solution), expected);
    check_matrix(assemble_spatial<mechanical::stiffness_assembler_1d<T, I>>(
                     mesh, mechanical_parameters, settings, solution), expected);
}

const suite<"assemblers_1d"> _ = [] {
    "local_spatial_operator"_test = [] {
        // c = 2 on [0, 2], with two linear elements of length 1.
        // Symmetric K: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected = { { { 2., -2.,  0. }, //
                                                 { 0.,  4., -2. }, //
                                                 { 0.,  0.,  2. } } };
        check_spatial(make_mesh(), T{ 2 }, {}, { .theories = { theory_t::LOCAL } }, Expected);
    };

    "constant_mass_and_capacity"_test = [] {
        const auto mesh = make_mesh();
        const thermal::parameters_1d<T> thermal_parameters{
            { .physical = { .conductivity = T{ 1 }, .capacity = T{ 2 }, .density = T{ 3 } } }
        };
        const mechanical::parameters_1d<T> mechanical_parameters{
            { .physical = { .youngs_modulus = T{ 1 }, .density = T{ 6 } } }
        };
        // The mass factor is rho*c_p = 3*2 = 6 for heat and rho = 6 for mechanics.
        // Symmetric M: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected = { { { 2., 1., 0. }, //
                                                 { 0., 4., 1. }, //
                                                 { 0., 0., 2. } } };
        check_matrix(assemble_mass<thermal::heat_capacity_assembler_1d<T, I>>(mesh, thermal_parameters, {}), Expected);
        check_matrix(assemble_mass<mechanical::mass_assembler_1d<T, I>>(mesh, mechanical_parameters, {}), Expected);
    };

    "composite_segments"_test = [] {
        const auto mesh = std::make_shared<mesh::mesh_1d<T>>(
            metamath::finite_element::make_element_1d_integrated<T>(1, 2),
            std::vector<mesh::segment_data<T>>{ { .length = T{ 1 }, .elements = 1 },
                                                { .length = T{ 1 }, .elements = 1 } });
        const thermal::parameters_1d<T> thermal_parameters{
            { .physical = { .conductivity = T{ 2 } } }, { .physical = { .conductivity = T{ 4 } } }
        };
        const mechanical::parameters_1d<T> mechanical_parameters{
            { .physical = { .youngs_modulus = T{ 2 }, .density = T{ 1 } } },
            { .physical = { .youngs_modulus = T{ 4 }, .density = T{ 1 } } }
        };
        const problem_settings settings{ .theories = { theory_t::LOCAL, theory_t::LOCAL } };
        // c = 2 on [0, 1] and c = 4 on [1, 2]; the middle diagonal is 2 + 4.
        // Symmetric K: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected = { { { 2., -2.,  0. }, //
                                                 { 0.,  6., -4. }, //
                                                 { 0.,  0.,  4. } } };
        check_matrix(assemble_spatial<thermal::thermal_conductivity_assembler_1d<T, I>>(
                         mesh, thermal_parameters, settings), Expected);
        check_matrix(assemble_spatial<mechanical::stiffness_assembler_1d<T, I>>(
                         mesh, mechanical_parameters, settings), Expected);
    };

    "variable_coefficients"_test = [] {
        const auto mesh = make_mesh();
        const coefficient_t<T, 1> coefficient = spatial_dependency<T, 1>{
            [](const point<T, 1>& x) constexpr noexcept { return T{ 1 } + x; }
        };
        // c(x) = 1 + x gives element stiffness factors 3/2 and 5/2.
        // Symmetric K: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected_Spatial = { { { 1.5, -1.5,   0. }, //
                                                         {  0.,   4., -2.5 }, //
                                                         {  0.,   0.,  2.5 } } };
        check_spatial(mesh, coefficient, {}, { .theories = { theory_t::LOCAL } }, Expected_Spatial);
        const mechanical::parameters_1d<T> parameters{
            { .physical = { .youngs_modulus = T{ 1 }, .density = coefficient } }
        };
        // m(x) = 1 + x; mass entries are integrals of m(x)*N_i(x)*N_j(x).
        // Symmetric M: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected_Mass = { { { 5. / 12., 1. / 4.,        0. }, //
                                                      {       0., 4. / 3.,  5. / 12. }, //
                                                      {       0.,      0., 11. / 12. } } };
        check_matrix(assemble_mass<mechanical::mass_assembler_1d<T, I>>(mesh, parameters, {}), Expected_Mass);
    };

    "solution_dependent_coefficient"_test = [] {
        const auto mesh = make_mesh();
        const coefficient_t<T, 1> coefficient = solution_dependency<T, 1>{
            [](const point<T, 1>& x, const T u) constexpr noexcept { return x + u; }
        };
        // u = 1 at all quadrature nodes, so c(x, u) = x + u = 1 + x.
        // Both K and M therefore match variable_coefficients above.
        // Both are symmetric: check only their upper triangles stored in inner.
        static constexpr matrix_3 Expected = { { { 1.5, -1.5,   0. }, //
                                                 {  0.,   4., -2.5 }, //
                                                 {  0.,   0.,  2.5 } } };
        const std::vector<T> solution(mesh->qnodes_count(), T{ 1 });
        check_spatial(mesh, coefficient, {}, { .theories = { theory_t::LOCAL }, .is_solution_dependent = true },
                        Expected, solution);
        const mechanical::parameters_1d<T> parameters{
            { .physical = { .youngs_modulus = T{ 1 }, .density = coefficient } }
        };
        finite_element_matrix_1d<T, I> mass;
        init_matrix_portrait(mass.inner, *mesh, { .theories = { theory_t::LOCAL } });
        mechanical::mass_assembler_1d<T, I>{ mass, mesh }.calc_matrix(parameters, {}, solution);
        static constexpr matrix_3 Expected_Mass = { { { 5. / 12., 1. / 4.,        0. }, //
                                                      {       0., 4. / 3.,  5. / 12. }, //
                                                      {       0.,      0., 11. / 12. } } };
        check_matrix(mass, Expected_Mass);
    };

    "nonlocal_spatial_operator"_test = [] {
        // phi(x, y) = 1 is intentionally unnormalized: its integral over [0, 2] is 2.
        // This simple kernel makes the nonlocal integrals easy to check by hand.
        const model_parameters<1, T> model{
            .influence = [](const T&, const T&) constexpr noexcept { return T{ 1 }; }, .local_weight = T{ 0.5 }
        };
        // With c = 2 and local weight 1/2: K = K_local/2 + d*d^T, where d = (-1, 0, 1).
        // Symmetric K: check only the upper triangle stored in inner.
        static constexpr matrix_3 Expected = { { { 2., -1., -1. }, //
                                                 { 0.,  2., -1. }, //
                                                 { 0.,  0.,  2. } } };
        check_spatial(make_mesh(), T{ 2 }, model, { .theories = { theory_t::NONLOCAL } }, Expected);
    };

    "nonsymmetric_nonlocal_operator"_test = [] {
        const model_parameters<1, T> model{
            .influence = [](const T&, const T&) constexpr noexcept { return T{ 1 }; }, .local_weight = T{ 0.5 }
        };
        const coefficient_t<T, 1> coefficient = spatial_dependency<T, 1>{
            [](const point<T, 1>& x) constexpr noexcept { return T{ 1 } + x; }
        };
        // Same unnormalized kernel; c(y) = 1 + y multiplies the column basis derivative.
        // K is nonsymmetric: check the full matrix stored in inner.
        static constexpr matrix_3 Expected = { { {   1.5, -0.25, -1.25 }, //
                                                 { -0.75,    2., -1.25 }, //
                                                 { -0.75, -1.75,   2.5 } } };
        check_spatial(make_mesh(), coefficient, model,
                        { .theories = { theory_t::NONLOCAL }, .is_nonconstant_parameters = true }, Expected);
    };

    "neumann_integral_constraint"_test = [] {
        const auto mesh = make_mesh();
        const thermal::parameters_1d<T> parameters{ { .physical = { .conductivity = T{ 2 } } } };
        // Neumann constraint: integral(u) = Q; b_i = integral(N_i) = (1/2, 1, 1/2).
        // The full symmetric matrix is [K b; b^T 0]; check its upper triangle stored in inner.
        static constexpr std::array<std::array<T, 4>, 4> Expected = { { { 2., -2.,  0., 0.5 }, //
                                                                        { 0.,  4., -2.,  1. }, //
                                                                        { 0.,  0.,  2., 0.5 }, //
                                                                        { 0.,  0.,  0.,  0. } } };
        check_matrix(assemble_spatial<thermal::thermal_conductivity_assembler_1d<T, I>>(
                         mesh, parameters, { .theories = { theory_t::LOCAL }, .is_neumann = true }), Expected);
    };

    "dirichlet_after_matrix_combination"_test = [] {
        const auto mesh = make_mesh();
        const mechanical::parameters_1d<T> parameters{
            { .physical = { .youngs_modulus = T{ 2 }, .density = T{ 6 } } }
        };
        static constexpr std::array<bool, 2> Is_First_Kind = { true, true };
        const problem_settings settings{ .theories = { theory_t::LOCAL }, .is_first_kind = Is_First_Kind };
        auto stiffness = assemble_spatial<mechanical::stiffness_assembler_1d<T, I>>(mesh, parameters, settings);
        const auto mass = assemble_mass<mechanical::mass_assembler_1d<T, I>>(mesh, parameters, Is_First_Kind);
        // Dirichlet conditions at nodes 0 and 2; their diagonals stay zero until K - 4*M is formed.
        // K and M are symmetric: check their upper triangles stored in inner.
        static constexpr matrix_3 Expected_Stiffness = { { { 0., 0., 0. }, //
                                                           { 0., 4., 0. }, //
                                                           { 0., 0., 0. } } };
        static constexpr matrix_3 Expected_Mass = Expected_Stiffness;
        check_matrix(stiffness, Expected_Stiffness);
        check_matrix(mass, Expected_Mass);
        // Couplings between node 1 and each boundary are stored separately in bound.
        for (const size_t boundary : std::ranges::iota_view{ 0zu, 2zu }) {
            expect(eq(stiffness.bound[boundary].size(), 1zu)) << fatal;
            expect(eq(mass.bound[boundary].size(), 1zu)) << fatal;
            expect(approx(stiffness.bound[boundary].at(1), T{ -2 }, T{ 1e-14 }));
            expect(approx(mass.bound[boundary].at(1), T{ 1 }, T{ 1e-14 }));
        }
        stiffness.inner -= T{ 4 } * mass.inner;
        stiffness.set_first_kind_diagonal(Is_First_Kind);
        // Set the Dirichlet diagonals to 1 after combination; the middle entry is 4 - 4*4 = -12.
        // The result is symmetric: check its upper triangle stored in inner.
        static constexpr matrix_3 Expected = { { { 1.,   0., 0. }, //
                                                 { 0., -12., 0. }, //
                                                 { 0.,   0., 1. } } };
        check_matrix(stiffness, Expected);
    };
};

} // namespace

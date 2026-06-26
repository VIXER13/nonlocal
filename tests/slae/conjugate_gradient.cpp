#include <solvers/slae_metamath_powered/conjugate_gradient.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

suite<"conjugate_gradient"> _ = [] {
    "conjugate_gradient"_test = [] {
        // Poisson equation in 1D with Dirichlet boundary conditions
        // (1  0  0  0  0  0  0  0  0  0)
        // (0  4 -1  0  0  0  0  0  0  0)
        // (0  0  4 -1  0  0  0  0  0  0)
        // (0  0  0  4 -1  0  0  0  0  0)
        // (0  0  0  0  4 -1  0  0  0  0)
        // (0  0  0  0  0  4 -1  0  0  0)
        // (0  0  0  0  0  0  4 -1  0  0)
        // (0  0  0  0  0  0  0  4 -1  0)
        // (0  0  0  0  0  0  0  0  4  0)
        // (0  0  0  0  0  0  0  0  0  1)
        sparse_matrix<T> matrix{10, 10};
        matrix.portrait.shifts = {0, 1, 3, 5, 7, 9, 11, 13, 15, 16, 17};
        matrix.portrait.indices = {0, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9};
        matrix.values = {1., 4., -1., 4., -1., 4., -1., 4., -1., 4., -1., 4., -1., 4., -1., 4., 1.};

        const std::vector<T> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        const std::vector<T> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const auto solver = conjugate_gradient{matrix};
        const std::vector<T> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Conjugate gradient solver failed with diff = " << diff;
    };

    "conjugate_gradient_block"_test = [] {
        // Poisson equation in 1D with Dirichlet boundary conditions in blocks representation
        // ((1  0) | ( 0  0) | ( 0  0) | ( 0  0) | ( 0  0))
        // ((0  4) | (-1  0) | ( 0  0) | ( 0  0) | ( 0  0))
        // ------------------------------------------------
        // ((0  0) | ( 4 -1) | ( 0  0) | ( 0  0) | ( 0  0))
        // ((0  0) | (-1  4) | (-1  0) | ( 0  0) | ( 0  0))
        // ------------------------------------------------
        // ((0  0) | ( 0  0) | ( 4 -1) | ( 0  0) | ( 0  0))
        // ((0  0) | ( 0  0) | (-1  4) | (-1  0) | ( 0  0))
        // ------------------------------------------------
        // ((0  0) | ( 0  0) | ( 0  0) | ( 4 -1) | ( 0  0))
        // ((0  0) | ( 0  0) | ( 0  0) | (-1  4) | (-1  0))
        // ------------------------------------------------
        // ((0  0) | ( 0  0) | ( 0  0) | ( 0  0) | ( 4  0))
        // ((0  0) | ( 0  0) | ( 0  0) | ( 0  0) | ( 0  1))
        sparse_matrix<square_matrix<T, 2>> block_matrix{5, 5};
        block_matrix.portrait.shifts = {0, 2, 4, 6, 8, 9};
        block_matrix.portrait.indices = {0, 1, 1, 2, 2, 3, 3, 4, 4};
        block_matrix.values = {{1.,  0.,  0., 4.}, {0., 0., -1., 0.},
                               {4., -1., -1., 4.}, {0., 0., -1., 0.},
                               {4., -1., -1., 4.}, {0., 0., -1., 0.},
                               {4., -1., -1., 4.}, {0., 0., -1., 0.},
                               {4.,  0.,  0., 1.}};

        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = block_matrix.self_adjoint<matrix_part::Upper>() * expected;
        const auto solver = conjugate_gradient{block_matrix};
        const std::vector<std::array<T, 2>> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Conjugate gradient solver failed with diff = " << diff;
    };
};

}
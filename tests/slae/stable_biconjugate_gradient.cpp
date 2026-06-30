#include <solvers/slae_metamath_powered/stable_biconjugate_gradient.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

suite<"stable_biconjugate_gradient"> _ = [] {
    "scalar"_test = [] {
        // Asymmetrical matrix for testing:
        // [1  0  0  0  0  0  0  0  0  0]
        // [0  4 -1  0  0  0  0  0  0  0]
        // [0  2  4  1  0  0  0  0  0  0]
        // [0  0 -2  4 -1  0  0  0  0  0]
        // [0  0  0 -1  4 -1  0  0  0  0]
        // [0  0  0  0 -2  4 -1  0  0  0]
        // [0  0  0  0  0  0  4 -1  0  0]
        // [0  0  0  0  0  0  2  4 -1  0]
        // [0  0  0  0  0  0  0  2  4  0]
        // [0  0  0  0  0  0  0  0  0  1]
        sparse_matrix<T> matrix{10, 10};
        matrix.portrait.shifts = {0, 1, 3, 6, 9, 12, 15, 17, 20, 22, 23};
        matrix.portrait.indices = {0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 6, 7, 6, 7, 8, 7, 8, 9};
        matrix.values = {1., 4., -1., 2., 4., 1., -2., 4., -1., -1., 4., -1., -2., 4., -1., 4., -1., 2., 4., -1., 2., 4., 1.};

        const std::vector<T> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        const std::vector<T> b = matrix * expected;
        const auto solver = stable_biconjugate_gradient{matrix};
        const std::vector<T> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };

    "block"_test = [] {
        // Asymmetrical matrix for testing in blocks representation:
        // [ 1  0 |  0  0 |  0  0 |  0  0 |  0  0 ]
        // [ 0  4 | -1  0 |  0  0 |  0  0 |  0  0 ]
        // ----------------------------------------
        // [ 0  2 |  4  1 |  0  0 |  0  0 |  0  0 ]
        // [ 0  0 | -2  4 | -1  0 |  0  0 |  0  0 ]
        // ----------------------------------------
        // [ 0  0 |  0 -1 |  4 -1 |  0  0 |  0  0 ]
        // [ 0  0 |  0  0 | -2  4 | -1  0 |  0  0 ]
        // ----------------------------------------
        // [ 0  0 |  0  0 |  0  0 |  4 -1 |  0  0 ]
        // [ 0  0 |  0  0 |  0  0 |  2  4 | -1  0 ]
        // ----------------------------------------
        // [ 0  0 |  0  0 |  0  0 |  0  2 |  4  0 ]
        // [ 0  0 |  0  0 |  0  0 |  0  0 |  0  1 ]
        sparse_matrix<square_matrix<T, 2>> block_matrix{5,5};
        block_matrix.portrait.shifts = {0, 2, 5, 8, 10, 12};
        block_matrix.portrait.indices = {0, 1, 0, 1, 2, 1, 2, 3, 3, 4, 3, 4};
        block_matrix.values = {{1.,  0., 0., 4.}, {0.,  0., -1., 0.},
                               {0.,  2., 0., 0.}, {4.,  1., -2., 4.}, {0., 0., -1., 0.},
                               {0., -1., 0., 0.}, {4., -1., -2., 4.}, {0., 0., -1., 0.},
                               {4., -1., 2., 4.}, {0., 0., -1., 0.},
                               {0.,  2., 0., 0.}, {4.,  0., 0., 1.}};

        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = block_matrix * expected;
        const auto solver = stable_biconjugate_gradient{block_matrix};
        const std::vector<std::array<T, 2>> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };
};

}
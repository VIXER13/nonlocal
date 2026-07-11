#include <metamath/linear/linear.hpp>
#include <solvers/slae_metamath_powered/diagonal_preconditioner.hpp>
#include <solvers/slae_metamath_powered/identity_preconditioner.hpp>
#include <solvers/slae_metamath_powered/ildlt_preconditioner.hpp>
#include <solvers/slae_metamath_powered/ilu0_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Epsilon = std::numeric_limits<T>::epsilon();
constexpr auto Inf = metamath::constants::Infinity<size_t>;

// 10x10 symmetric matrix
// [10 -1 -0.5    0    0    0    0    0    0    0]
// [ 0 10   -1 -0.5    0    0    0    0    0    0]
// [ 0  0   10   -1 -0.5    0    0    0    0    0]
// [ 0  0    0   10   -1 -0.5    0    0    0    0]
// [ 0  0    0    0   10   -1 -0.5    0    0    0]
// [ 0  0    0    0    0   10   -1 -0.5    0    0]
// [ 0  0    0    0    0    0   10   -1 -0.5    0]
// [ 0  0    0    0    0    0    0   10   -1 -0.5]
// [ 0  0    0    0    0    0    0    0   10   -1]
// [ 0  0    0    0    0    0    0    0    0   10]
sparse_matrix<T> scalar_symmetric_matrix() {
    sparse_matrix<T> matrix{10, 10};
    matrix.portrait.shifts  = {0, 3, 6, 9, 12, 15, 18, 21, 24, 26, 27};
    matrix.portrait.indices = {0, 1, 2,  1, 2, 3,  2, 3, 4,  3, 4, 5, 
                               4, 5, 6,  5, 6, 7,  6, 7, 8,  7, 8, 9, 
                                  8, 9,        9};
    matrix.values           = {10., -1., -0.5,  10., -1., -0.5,  10., -1., -0.5,  10., -1., -0.5,
                               10., -1., -0.5,  10., -1., -0.5,  10., -1., -0.5,  10., -1., -0.5,
                               10., -1.,        10.};
    return matrix;
}

// 10x10 general matrix
// [10 -1 -0.5    0    0    0    0    0    0    0]
// [ 2 10   -1 -0.5    0    0    0    0    0    0]
// [ 1  2   10   -1 -0.5    0    0    0    0    0]
// [ 0  1    2   10   -1 -0.5    0    0    0    0]
// [ 0  0    1    2   10   -1 -0.5    0    0    0]
// [ 0  0    0    1    2   10   -1 -0.5    0    0]
// [ 0  0    0    0    1    2   10   -1 -0.5    0]
// [ 0  0    0    0    0    1    2   10   -1 -0.5]
// [ 0  0    0    0    0    0    1    2   10   -1]
// [ 0  0    0    0    0    0    0    1    2   10]
sparse_matrix<T> scalar_general_matrix() {
    sparse_matrix<T> matrix{10, 10};
    matrix.portrait.shifts  = {0, 3, 7, 12, 17, 22, 27, 32, 37, 41, 44};
    matrix.portrait.indices = {
              0, 1, 2,     0, 1, 2, 3,
        0, 1, 2, 3, 4,  1, 2, 3, 4, 5,   2, 3, 4, 5, 6,
        3, 4, 5, 6, 7,  4, 5, 6, 7, 8,  5, 6, 7, 8, 9,
        6, 7, 8, 9,     7, 8, 9
    };
    matrix.values = {
                10., -1., -0.5,      2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,  1., 2., 10., -1., -0.5,  1., 2., 10., -1., -0.5,
        1., 2., 10., -1., -0.5,  1., 2., 10., -1., -0.5,  1., 2., 10., -1., -0.5,
        1., 2., 10., -1.,        1., 2., 10.
    };
    return matrix;
}

// 5x5 symmetric block matrix
// [ 10   -1 |  0    0 |  0    0 |  0    0 |  0    0 ]
// [ -1   10 | -1    0 |  0 -0.5 |  0    0 |  0    0 ]
// ---------------------------------------------------
// [  0    0 | 10   -1 |  0    0 |  0    0 |  0    0 ]
// [  0    0 | -1   10 | -1    0 |  0 -0.5 |  0    0 ]
// ---------------------------------------------------
// [  0    0 |  0    0 | 10   -1 |  0    0 |  0    0 ]
// [  0    0 |  0    0 | -1   10 | -1    0 |  0 -0.5 ]
// ---------------------------------------------------
// [  0    0 |  0    0 |  0    0 | 10   -1 |  0    0 ]
// [  0    0 |  0    0 |  0    0 | -1   10 | -1    0 ]
// ---------------------------------------------------
// [  0    0 |  0    0 |  0    0 |  0    0 | 10    0 ]
// [  0    0 |  0    0 |  0    0 |  0    0 |  0   10 ]
sparse_matrix<square_matrix<T, 2>> block_symmetric_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{5, 5};
    matrix.portrait.shifts  = {0, 3, 6, 9, 11, 12};
    matrix.portrait.indices = {0, 1, 2,  1, 2, 3,  2, 3, 4,  3, 4,  4};
    matrix.values = { {10., -1., -1., 10.}, {0.,0.,-1.,0.}, {0.,0.,0.,-.5},
                      {10., -1., -1., 10.}, {0.,0.,-1.,0.}, {0.,0.,0.,-.5},
                      {10., -1., -1., 10.}, {0.,0.,-1.,0.}, {0.,0.,0.,-.5},
                      {10., -1., -1., 10.}, {0.,0.,-1.,0.},
                      {10., -1., -1., 10.} };
    return matrix;
}

// 5x5 general block matrix
// [ 10   -1 |   0    0 |   0    0 |  0    0 |  0    0 ]
// [ -1   10 |  -1    0 |   0 -0.5 |  0    0 |  0    0 ]
// ---------------------------------------------------
// [  0    2 |  10   -1 |   0    0 |  0    0 |  0    0 ]
// [  0    0 |  -1   10 |  -1    0 |  0 -0.5 |  0    0 ]
// ---------------------------------------------------
// [  0    0 | 0.5    2 |  10   -1 |  0    0 |  0    0 ]
// [  0    0 |   0    0 |  -1   10 | -1    0 |  0 -0.5 ]
// ---------------------------------------------------
// [  0    0 | 0.5    0 |   0    2 | 10   -1 |  0    0 ]
// [  0    0 |   0    0 |   0    0 | -1   10 | -1    0 ]
// ---------------------------------------------------
// [  0    0 |   0    0 | 0.5    0 |  0    2 | 10    0 ]
// [  0    0 |   0    0 |   0    0 |  0    0 |  0   10 ]
sparse_matrix<square_matrix<T, 2>> block_general_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{5, 5};
    matrix.portrait.shifts  = {0, 3, 7, 12, 16, 19};
    matrix.portrait.indices = {0, 1, 2,  0, 1, 2, 3,  0, 1, 2, 3, 4,  1, 2, 3, 4,  2, 3, 4};
    matrix.values = {                        {10., -1., -1., 10.}, {0., 0., -1., 0.}, {0., 0., 0., -0.5},
                           {0., 2., 0., 0.}, {10., -1., -1., 10.}, {0., 0., -1., 0.}, {0., 0., 0., -0.5},
        {0.5, 0., 0., 0.}, {0., 2., 0., 0.}, {10., -1., -1., 10.}, {0., 0., -1., 0.}, {0., 0., 0., -0.5},
        {0.5, 0., 0., 0.}, {0., 2., 0., 0.}, {10., -1., -1., 10.}, {0., 0., -1., 0.},
        {0.5, 0., 0., 0.}, {0., 2., 0., 0.}, {10., -1., -1., 10.} };
    return matrix;
}

// -------------------------------------------------------------------

suite<"identity_preconditioner"> _identity = [] {
    "scalar"_test = [] {
        const identity_preconditioner<T> preconditioner;
        const std::vector<T> r = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        expect(preconditioner.solve(r) == r) << "identity preconditioner must return input unchanged";
    };

    "block"_test = [] {
        const identity_preconditioner<square_matrix<T, 2>> preconditioner;
        const std::vector<std::array<T, 2>> r = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        expect(preconditioner.solve(r) == r) << "block identity preconditioner must return input unchanged";
    };
};

// -------------------------------------------------------------------

suite<"diagonal_preconditioner"> _diagonal = [] {
    "scalar"_test = [] {
        // All diagonal entries are 10, so inv(D)*r = r/10.
        const diagonal_preconditioner<T> preconditioner{scalar_symmetric_matrix()};
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> r = 10 * expected;
        const T diff = norm<Inf>(preconditioner.solve(r) - expected);
        expect(approx(diff, 0.0, Epsilon)) << "diagonal scalar solve failed, diff=" << diff;
    };

    "block"_test = [] {
        // Diagonal blocks are all [[10,-1],[-1,10]].
        const diagonal_preconditioner<square_matrix<T, 2>> preconditioner{block_symmetric_matrix()};
        static constexpr square_matrix<T, 2> diag_block = {10., -1., -1., 10.};
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        std::vector<std::array<T, 2>> r(expected.size());
        for (const size_t i : std::ranges::iota_view{0zu, expected.size()})
            r[i] = diag_block * expected[i];
        const T diff = norm<Inf>(preconditioner.solve(r) - expected);
        expect(approx(diff, 0.0, Epsilon)) << "diagonal block solve failed, diff=" << diff;
    };
};

// -------------------------------------------------------------------

suite<"ildlt_preconditioner"> _ildlt = [] {
    "scalar_factorization"_test = [] {
        // ILU0/ILDLT is exact for banded matrices: (L D L^T)^{-1} A x = x.
        const auto matrix = scalar_symmetric_matrix();
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{scalar_symmetric_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 8.9e-16;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt scalar exact solve failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        const auto matrix = block_symmetric_matrix();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{block_symmetric_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt block exact solve failed, diff=" << diff;
    };
};

// -------------------------------------------------------------------

suite<"ilu0_preconditioner"> _ilu0 = [] {
    "scalar_factorization"_test = [] {
        const auto matrix = scalar_general_matrix();
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> b = matrix * expected;
        const ilu0_preconditioner preconditioner{scalar_general_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ilu0 scalar exact solve failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        const auto matrix = block_general_matrix();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = matrix * expected;
        const ilu0_preconditioner preconditioner{block_general_matrix()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ilu0 block exact solve failed, diff=" << diff;
    };
};

}
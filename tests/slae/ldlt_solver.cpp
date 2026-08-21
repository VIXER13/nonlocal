#include <metamath/linear/linear.hpp>
#include <solvers/slae/ldlt_solver.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

// 6x6 symmetric sparse matrix, upper triangle only.
// Connections: (0,3), (1,4), (2,5), (0,5) — no adjacent off-diagonals, requires fill-in.
// Full symmetric view:
// [20  0  0 -3  0 -1]
// [ 0 18  0  0 -2  0]
// [ 0  0 16  0  0 -4]
// [-3  0  0 20  0  0]
// [ 0 -2  0  0 18  0]
// [-1  0 -4  0  0 16]
sparse_matrix<T> ldlt_scalar_matrix() {
    sparse_matrix<T> matrix{6, 6};
    matrix.portrait.shifts  = {0, 3, 5, 7, 8, 9, 10};
    matrix.portrait.indices = {0, 3, 5,  1, 4,  2, 5,  3,  4,  5};
    matrix.values           = {20., -3., -1.,  18., -2.,  16., -4.,  20.,  18.,  16.};
    return matrix;
}

// 4x4 symmetric block matrix (2x2 blocks), upper triangle only.
// Block connections: (0,2), (1,3) only — diagonal blocks plus one skip-1 off-diagonal pair.
// No adjacent block coupling at all, requires fill-in in block LDLT.
sparse_matrix<square_matrix<T, 2>> ldlt_block_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{4, 4};
    matrix.portrait.shifts  = {0, 2, 3, 5, 6};
    matrix.portrait.indices = {0, 2,  1,  2, 3,  3};
    // Block (0,0)=20I+J, (0,2)=-2I, (1,1)=18I+J, (2,2)=20I+J, (2,3)=-2I, (3,3)=18I+J
    matrix.values = {
        {20., 1., 1., 20.}, {-2., 0., 0., -2.},
        {18., 1., 1., 18.},
        {20., 1., 1., 20.}, {-2., 0., 0., -2.},
        {18., 1., 1., 18.}
    };
    return matrix;
}

suite<"ldlt_solver"> _ldlt = [] {
    "scalar_factorization"_test = [] {
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6.};
        const std::vector<T> b = ldlt_scalar_matrix().self_adjoint<matrix_part::Upper>() * expected;
        const ldlt_solver solver{ldlt_scalar_matrix()};
        const T diff = norm<Inf>(solver.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-13;
        expect(approx(diff, 0.0, Epsilon)) << "ldlt scalar exact solve failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}};
        const std::vector<std::array<T, 2>> b = ldlt_block_matrix().self_adjoint<matrix_part::Upper>() * expected;
        const ldlt_solver solver{ldlt_block_matrix()};
        const T diff = norm<Inf>(solver.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-13;
        expect(approx(diff, 0.0, Epsilon)) << "ldlt block exact solve failed, diff=" << diff;
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { ldlt_solver{sparse_matrix<T>{3, 4}}; })) <<
            "ldlt solver must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const ldlt_solver solver{ldlt_scalar_matrix()};
        const std::vector<T> rhs{1., 2.};
        expect(throws<std::invalid_argument>([&solver, &rhs] { solver.solve(rhs); })) <<
            "ldlt solver must throw for rhs vector of wrong size";
    };
};

}

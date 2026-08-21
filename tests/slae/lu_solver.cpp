#include <metamath/linear/linear.hpp>
#include <solvers/slae/lu_solver.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

// 6x6 general sparse matrix.
// Non-zeros: diagonal + (0,3), (0,4), (2,5), (3,0), (5,1), (5,2) — scattered pattern, requires fill-in.
// [ 15  0  0 -2 -1  0]
// [  0 14  0  0  0  0]
// [  0  0 12  0  0 -3]
// [ -4  0  0 15  0  0]
// [  0  0  0  0 14  0]
// [  0  2  1  0  0 12]
sparse_matrix<T> lu_scalar_matrix() {
    sparse_matrix<T> matrix{6, 6};
    matrix.portrait.shifts  = {0, 3, 4, 6, 8, 9, 12};
    matrix.portrait.indices = {0, 3, 4,  1,  2, 5,  0, 3,  4,  1, 2, 5};
    matrix.values           = {15., -2., -1.,  14.,  12., -3.,  -4., 15.,  14.,  2., 1., 12.};
    return matrix;
}

// 4x4 general block matrix (2x2 blocks).
// Block non-zeros: diagonal + (0,2), (1,3), (2,0), (3,1) — no adjacent coupling, requires fill-in.
sparse_matrix<square_matrix<T, 2>> lu_block_matrix() {
    sparse_matrix<square_matrix<T, 2>> matrix{4, 4};
    matrix.portrait.shifts  = {0, 2, 3, 5, 6};
    matrix.portrait.indices = {0, 2,  1,  0, 2,  3};
    matrix.values = {
        {16., 1., 0., 16.}, {-2., 0., 0., -2.},
        {14., 1., 0., 14.},
        {3., 0., 0., 3.}, {16., 1., 0., 16.},
        {14., 1., 0., 14.}
    };
    return matrix;
}

suite<"lu_solver"> _lu = [] {
    "scalar_factorization"_test = [] {
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6.};
        const std::vector<T> b = lu_scalar_matrix() * expected;
        const lu_solver solver{lu_scalar_matrix()};
        const T diff = norm<Inf>(solver.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-13;
        expect(approx(diff, 0.0, Epsilon)) << "lu scalar exact solve failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}};
        const std::vector<std::array<T, 2>> b = lu_block_matrix() * expected;
        const lu_solver solver{lu_block_matrix()};
        const T diff = norm<Inf>(solver.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-13;
        expect(approx(diff, 0.0, Epsilon)) << "lu block exact solve failed, diff=" << diff;
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { lu_solver{sparse_matrix<T>{5, 3}}; })) <<
            "lu solver must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const lu_solver solver{lu_scalar_matrix()};
        const std::vector<T> rhs{1., 2., 3.};
        expect(throws<std::invalid_argument>([&solver, &rhs] { solver.solve(rhs); })) <<
            "lu solver must throw for rhs vector of wrong size";
    };
};

}

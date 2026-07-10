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

constexpr auto Inf = metamath::constants::Infinity<size_t>;

// 5x5 symmetric tridiagonal matrix (upper triangular storage):
// [ 4 -1  0  0  0 ]
// [ 0  4 -1  0  0 ]
// [ 0  0  4 -1  0 ]
// [ 0  0  0  4 -1 ]
// [ 0  0  0  0  4 ]
sparse_matrix<T> make_scalar_symmetric() {
    sparse_matrix<T> m{5, 5};
    m.portrait.shifts  = {0, 2, 4, 6, 8, 9};
    m.portrait.indices = {0, 1, 1, 2, 2, 3, 3, 4, 4};
    m.values           = {4., -1., 4., -1., 4., -1., 4., -1., 4.};
    return m;
}

// 5x5 general tridiagonal matrix (full CSR):
// [ 4 -1  0  0  0 ]
// [ 2  4 -1  0  0 ]
// [ 0  2  4 -1  0 ]
// [ 0  0  2  4 -1 ]
// [ 0  0  0  2  4 ]
sparse_matrix<T> make_scalar_general() {
    sparse_matrix<T> m{5, 5};
    m.portrait.shifts  = {0, 2, 5, 8, 11, 13};
    m.portrait.indices = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
    m.values           = {4., -1., 2., 4., -1., 2., 4., -1., 2., 4., -1., 2., 4.};
    return m;
}

// Block 3x3 symmetric matrix (2x2 blocks, upper triangular storage):
// Scalar equivalent:
// [ 4 -1 |  0  0 |  0  0 ]
// [-1  4 | -1  0 |  0  0 ]
// [ 0  0 |  4 -1 | -1  0 ]
// [ 0  0 | -1  4 |  0 -1 ]
// [ 0  0 |  0  0 |  4  0 ]
// [ 0  0 |  0  0 |  0  4 ]
sparse_matrix<square_matrix<T, 2>> make_block_symmetric() {
    sparse_matrix<square_matrix<T, 2>> m{3, 3};
    m.portrait.shifts  = {0, 2, 4, 5};
    m.portrait.indices = {0, 1, 1, 2, 2};
    m.values = {
        {4., -1., -1., 4.}, {0., 0., -1., 0.},
        {4., -1., -1., 4.}, {0., 0., -1., 0.},
        {4.,  0.,  0., 4.}
    };
    return m;
}

// Block 3x3 general matrix (2x2 blocks, full CSR):
// Scalar equivalent:
// [ 4 -1 |  0  0 |  0  0 ]
// [ 2  4 | -1  0 |  0  0 ]
// [ 0  0 |  4 -1 | -1  0 ]
// [ 0  2 | -1  4 |  0 -1 ]
// [ 0  0 |  0  0 |  4  0 ]
// [ 0  0 |  0  2 |  0  4 ]
sparse_matrix<square_matrix<T, 2>> make_block_general() {
    sparse_matrix<square_matrix<T, 2>> m{3, 3};
    m.portrait.shifts  = {0, 2, 5, 7};
    m.portrait.indices = {0, 1, 0, 1, 2, 1, 2};
    m.values = {
        {4., -1., 2., 4.}, {0., 0., -1., 0.},
        {0., 0., 0., 2.},  {4., -1., -1., 4.}, {0., 0., -1., 0.},
        {0., 2., 0., 0.},  {4., 0., 0., 4.}
    };
    return m;
}

// -------------------------------------------------------------------

suite<"identity_preconditioner"> _identity = [] {
    "scalar"_test = [] {
        identity_preconditioner<T> p;
        const std::vector<T> r = {1., 2., 3., 4., 5.};
        expect(p.solve(r) == r) << "identity preconditioner must return input unchanged";
    };

    "block"_test = [] {
        identity_preconditioner<square_matrix<T, 2>> p;
        const std::vector<std::array<T, 2>> r = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        expect(p.solve(r) == r) << "block identity preconditioner must return input unchanged";
    };
};

// -------------------------------------------------------------------

suite<"diagonal_preconditioner"> _diagonal = [] {
    "scalar"_test = [] {
        // Diagonal of make_scalar_symmetric: [4, 4, 4, 4, 4]
        diagonal_preconditioner<T> p{make_scalar_symmetric()};
        const std::vector<T> r = {4., 8., 12., 16., 20.};
        const std::vector<T> x = p.solve(r);
        const std::vector<T> expected = {1., 2., 3., 4., 5.};
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-15)) << "diagonal preconditioner scalar solve failed, diff=" << diff;
    };

    "block"_test = [] {
        // Diagonal blocks of make_block_symmetric: [{4,-1,-1,4}, {4,-1,-1,4}, {4,0,0,4}]
        diagonal_preconditioner<square_matrix<T, 2>> p{make_block_symmetric()};
        // rhs = diag_block * x for known x
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        const square_matrix<T, 2> d0 = {4., -1., -1., 4.};
        const square_matrix<T, 2> d2 = {4.,  0.,  0., 4.};
        const std::vector<std::array<T, 2>> r = {d0 * expected[0], d0 * expected[1], d2 * expected[2]};
        const std::vector<std::array<T, 2>> x = p.solve(r);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-14)) << "diagonal preconditioner block solve failed, diff=" << diff;
    };
};

// -------------------------------------------------------------------

suite<"ildlt_preconditioner"> _ildlt = [] {
    "scalar_factorization"_test = [] {
        // For a tridiagonal symmetric positive-definite matrix the ILDLT
        // factorization is exact (no fill-in exists) so (L D L^T)^{-1} b = A^{-1} b.
        auto p = ildlt_preconditioner{make_scalar_symmetric()};
        const std::vector<T> expected = {1., 2., 3., 4., 5.};
        const std::vector<T> b = make_scalar_symmetric().self_adjoint<matrix_part::Upper>() * expected;
        const std::vector<T> x = p.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-13)) << "ildlt scalar exact solve failed, diff=" << diff;
    };

    "scalar_symmetry"_test = [] {
        // For a symmetric problem and exact factorization: L D L^T x = b => x = A^{-1} b.
        // Verify the solution is self-consistent: A*(ILDLT^{-1} b) ≈ b.
        auto p = ildlt_preconditioner{make_scalar_symmetric()};
        const auto m = make_scalar_symmetric();
        const std::vector<T> b = {1., 2., 3., 4., 5.};
        const std::vector<T> x = p.solve(b);
        const std::vector<T> Ax = m.self_adjoint<matrix_part::Upper>() * x;
        const T diff = norm<Inf>(Ax - b);
        expect(approx(diff, 0.0, 1e-13)) << "ildlt scalar residual check failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        auto p = ildlt_preconditioner{make_block_symmetric()};
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        const std::vector<std::array<T, 2>> b = make_block_symmetric().self_adjoint<matrix_part::Upper>() * expected;
        const std::vector<std::array<T, 2>> x = p.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-13)) << "ildlt block exact solve failed, diff=" << diff;
    };

    "block_symmetry"_test = [] {
        auto p = ildlt_preconditioner{make_block_symmetric()};
        const auto m = make_block_symmetric();
        const std::vector<std::array<T, 2>> b = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        const std::vector<std::array<T, 2>> x = p.solve(b);
        const std::vector<std::array<T, 2>> Ax = m.self_adjoint<matrix_part::Upper>() * x;
        const T diff = norm<Inf>(Ax - b);
        expect(approx(diff, 0.0, 1e-13)) << "ildlt block residual check failed, diff=" << diff;
    };
};

// -------------------------------------------------------------------

suite<"ilu0_preconditioner"> _ilu0 = [] {
    "scalar_factorization"_test = [] {
        // For tridiagonal general matrix the ILU0 factorization is exact.
        auto p = ilu0_preconditioner{make_scalar_general()};
        const std::vector<T> expected = {1., 2., 3., 4., 5.};
        const std::vector<T> b = make_scalar_general() * expected;
        const std::vector<T> x = p.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-13)) << "ilu0 scalar exact solve failed, diff=" << diff;
    };

    "scalar_residual"_test = [] {
        auto p = ilu0_preconditioner{make_scalar_general()};
        const auto m = make_scalar_general();
        const std::vector<T> b = {1., 2., 3., 4., 5.};
        const std::vector<T> x = p.solve(b);
        const std::vector<T> Ax = m * x;
        const T diff = norm<Inf>(Ax - b);
        expect(approx(diff, 0.0, 1e-13)) << "ilu0 scalar residual check failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        auto p = ilu0_preconditioner{make_block_general()};
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        const std::vector<std::array<T, 2>> b = make_block_general() * expected;
        const std::vector<std::array<T, 2>> x = p.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1e-13)) << "ilu0 block exact solve failed, diff=" << diff;
    };

    "block_residual"_test = [] {
        auto p = ilu0_preconditioner{make_block_general()};
        const auto m = make_block_general();
        const std::vector<std::array<T, 2>> b = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}};
        const std::vector<std::array<T, 2>> x = p.solve(b);
        const std::vector<std::array<T, 2>> Ax = m * x;
        const T diff = norm<Inf>(Ax - b);
        expect(approx(diff, 0.0, 1e-13)) << "ilu0 block residual check failed, diff=" << diff;
    };
};

}

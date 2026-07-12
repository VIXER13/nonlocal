#include "create_matrix.hpp"

#include <solvers/slae_metamath_powered/ildlt_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;
    
constexpr auto Inf = metamath::constants::Infinity<size_t>;

suite<"ildlt_preconditioner"> _ildlt = [] {
    "scalar_factorization"_test = [] {
        // ILU0/ILDLT is exact for banded matrices: (L D L^T)^{-1} A x = x.
        const auto matrix = scalar_symmetric_matrix<T>();
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 8.9e-16;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt scalar exact solve failed, diff=" << diff;
    };

    "block_factorization"_test = [] {
        const auto matrix = block_symmetric_matrix<T>();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = matrix.self_adjoint<matrix_part::Upper>() * expected;
        const ildlt_preconditioner preconditioner{block_symmetric_matrix<T>()};
        const T diff = norm<Inf>(preconditioner.solve(b) - expected);
        static constexpr auto Epsilon = 1.8e-15;
        expect(approx(diff, 0.0, Epsilon)) << "ildlt block exact solve failed, diff=" << diff;
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { ildlt_preconditioner{sparse_matrix<T>{3, 4}}; })) <<
            "ildlt preconditioner must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const ildlt_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const std::vector<T> rhs{1., 2.};
        expect(throws<std::invalid_argument>([&preconditioner, &rhs] { preconditioner.solve(rhs); })) <<
            "ildlt preconditioner must throw for rhs vector of wrong size";
    };
};

}
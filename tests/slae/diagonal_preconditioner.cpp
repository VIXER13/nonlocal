#include "create_matrix.hpp"

#include <solvers/slae/diagonal_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Epsilon = std::numeric_limits<T>::epsilon();
constexpr auto Inf = metamath::constants::Infinity<size_t>;

suite<"diagonal_preconditioner"> _diagonal = [] {
    "scalar"_test = [] {
        // All diagonal entries are 10, so inv(D)*r = r/10.
        const diagonal_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const std::vector<T> expected = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10.};
        const std::vector<T> r = 10 * expected;
        const T diff = norm<Inf>(preconditioner.solve(r) - expected);
        expect(approx(diff, 0.0, Epsilon)) << "diagonal scalar solve failed, diff=" << diff;
    };

    "block"_test = [] {
        // Diagonal blocks are all [[10,-1],[-1,10]].
        const diagonal_preconditioner preconditioner{block_symmetric_matrix<T>()};
        static constexpr square_matrix<T, 2> diag_block = {10., -1., -1., 10.};
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        std::vector<std::array<T, 2>> r(expected.size());
        for (const size_t i : std::ranges::iota_view{0zu, expected.size()})
            r[i] = diag_block * expected[i];
        const T diff = norm<Inf>(preconditioner.solve(r) - expected);
        expect(approx(diff, 0.0, Epsilon)) << "diagonal block solve failed, diff=" << diff;
    };

    "wrong_matrix_size"_test = [] {
        expect(throws<std::invalid_argument>([] { diagonal_preconditioner{sparse_matrix<T>{3, 4}}; })) <<
            "diagonal preconditioner must throw for non-square matrix";
    };

    "wrong_vector_size"_test = [] {
        const diagonal_preconditioner preconditioner{scalar_symmetric_matrix<T>()};
        const std::vector<T> rhs{1., 2.};
        expect(throws<std::invalid_argument>([&preconditioner, &rhs] { preconditioner.solve(rhs); })) <<
            "diagonal preconditioner must throw for rhs vector of wrong size";
    };
};

}
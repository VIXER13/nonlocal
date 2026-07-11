#include "create_matrix.hpp"

#include <solvers/slae_metamath_powered/ilu0_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

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
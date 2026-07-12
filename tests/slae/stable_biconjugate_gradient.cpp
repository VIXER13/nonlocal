#include "create_matrix.hpp"

#include <solvers/slae_metamath_powered/stable_biconjugate_gradient.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace nonlocal::tests;
using namespace metamath::linear;
using namespace metamath::operators;
using T = double;

constexpr auto Inf = metamath::constants::Infinity<size_t>;

suite<"stable_biconjugate_gradient"> _ = [] {
    "scalar"_test = [] {
        const auto matrix = scalar_general_matrix<T>();
        const std::vector<T> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        const std::vector<T> b = matrix * expected;
        const auto solver = stable_biconjugate_gradient{matrix};
        const std::vector<T> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };

    "block"_test = [] {
        const auto block_matrix = block_general_matrix<T>();
        const std::vector<std::array<T, 2>> expected = {{{1., 2.}}, {{3., 4.}}, {{5., 6.}}, {{7., 8.}}, {{9., 10.}}};
        const std::vector<std::array<T, 2>> b = block_matrix * expected;
        const auto solver = stable_biconjugate_gradient{block_matrix};
        const std::vector<std::array<T, 2>> x = solver.solve(b);
        const T diff = norm<Inf>(x - expected);
        expect(approx(diff, 0.0, 1.8e-15)) << "Stable BiConjugate gradient solver failed with diff = " << diff;
    };
};

}
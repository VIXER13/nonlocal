#include <metamath/linear/fixed_matrix.hpp>
#include <solvers/slae_metamath_powered/identity_preconditioner.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace nonlocal::slae;
using namespace metamath::linear;
using T = double;

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

}
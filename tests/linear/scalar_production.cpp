#include <metamath/linear/scalar_production.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;
using T = double;

constexpr auto Epsilon = std::numeric_limits<T>::epsilon();

const suite<"scalar_production"> _ = [] {
    "scalar_production"_test = []<class Container> {
        const Container x{1.0, 2.0, 3.0};
        const Container y{4.0, 5.0, 6.0};
        expect(approx(scalar_production(x, y), 32.0, Epsilon));
    } | std::tuple<std::array<T, 3>, std::vector<T>>{};

    "scalar_production_block"_test = []<class Container> {
        const Container x{{{1.0, 2.0}, {3.0, 4.0}}};
        const Container y{{{5.0, 6.0}, {7.0, 8.0}}};
        expect(approx(scalar_production(x, y), 70.0, Epsilon));
    } | std::tuple<std::array<std::array<T, 2>, 2>, std::vector<std::array<T, 2>>>{};
};

}
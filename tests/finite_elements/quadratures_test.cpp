#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

std::array<std::unique_ptr<quadrature_1d_base<T>>, 5> init_quadratures() {
    return {
        std::make_unique<quadrature_1d<T, gauss, 1>>(),
        std::make_unique<quadrature_1d<T, gauss, 2>>(),
        std::make_unique<quadrature_1d<T, gauss, 3>>(),
        std::make_unique<quadrature_1d<T, gauss, 4>>(),
        std::make_unique<quadrature_1d<T, gauss, 5>>()
    };
}

const suite<"quadrature_1d"> _ = [] {
    size_t order = 1;
    for(const auto& quadrature : init_quadratures()) {
        const std::string suffix = "_order_" + std::to_string(order);

        test("nodes_count" + suffix) = [&quadrature, order] {
            expect(eq(quadrature->nodes_count(), order)) << "Unexpected nodes count.";
        };

        test("boundaries" + suffix) = [&quadrature] {
            expect(eq(quadrature->boundary(side_1d::LEFT), T{-1})) << "Unexpected left boundary position.";
            expect(eq(quadrature->boundary(side_1d::RIGHT), T{1})) << "Unexpected right boundary position.";
        };

        test("weights_sum" + suffix) = [&quadrature] {
            const auto nodes = quadrature->nodes();
            const auto weight_summator = [&quadrature](const T sum, const size_t node) {
                return sum + quadrature->weight(node);
            };
            const T weights_sum = std::accumulate(nodes.begin(), nodes.end(), T{0}, weight_summator);
            const T length = quadrature->boundary(side_1d::RIGHT) - quadrature->boundary(side_1d::LEFT);
            expect(lt(std::abs(weights_sum - length), std::numeric_limits<T>::epsilon())) << "Unexpected weights sum.";
        };

        ++order;
    }
};

}
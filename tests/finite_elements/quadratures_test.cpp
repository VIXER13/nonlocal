#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

const suite<"quadrature_1d"> _ = [] {
    size_t order = 1;
    for(const size_t order : std::ranges::iota_view{1zu, 6zu}) {
        const std::string suffix = "_order_" + std::to_string(order);
        const auto quadrature = make_quadrature_1d<T>(order);

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

        test("copy" + suffix) = [&quadrature] {
            const auto copied_quadrature = quadrature->copy();
            expect(eq(copied_quadrature->nodes_count(), quadrature->nodes_count())) << "Unexpected nodes count in copied quadrature.";
            expect(eq(copied_quadrature->boundary(side_1d::LEFT), quadrature->boundary(side_1d::LEFT))) << "Unexpected left boundary position in copied quadrature.";
            expect(eq(copied_quadrature->boundary(side_1d::RIGHT), quadrature->boundary(side_1d::RIGHT))) << "Unexpected right boundary position in copied quadrature.";
            for(const size_t i : quadrature->nodes()) {
                expect(eq(copied_quadrature->weight(i), quadrature->weight(i))) << "Unexpected weight in copied quadrature.";
                expect(eq(copied_quadrature->node(i), quadrature->node(i))) << "Unexpected node in copied quadrature.";
            }
        };
    }
};

}
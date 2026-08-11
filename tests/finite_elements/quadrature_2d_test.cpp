#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

const suite<"quadrature_2d"> _ = [] {
    for(const size_t order : std::ranges::iota_view{1zu, 6zu}) {
        const std::string suffix = "_order_" + std::to_string(order);
        const auto quadrature = make_quadrature_2d<T>(order);

        test("nodes_count" + suffix) = [&quadrature, order] {
            expect(eq(quadrature->nodes_count(), order * order)) << "Unexpected nodes count.";
        };

        test("weights_sum" + suffix) = [&quadrature] {
            const auto nodes = quadrature->nodes();
            const auto weight_summator = [&quadrature](const T sum, const size_t node) {
                return sum + quadrature->weight(node);
            };
            const T weights_sum = std::accumulate(nodes.begin(), nodes.end(), T{0}, weight_summator);
            static constexpr T Epsilon = 9e-16;
            expect(approx(weights_sum, dynamic_cast<const geometry_2d_base<T>&>(*quadrature).area(), Epsilon)) << "Unexpected weights sum.";
        };

        test("copy" + suffix) = [&quadrature] {
            const auto copied_quadrature = quadrature->copy();
            expect(eq(copied_quadrature->nodes_count(), quadrature->nodes_count())) << "Unexpected nodes count in copied quadrature.";
            for(const size_t i : quadrature->nodes()) {
                expect(eq(copied_quadrature->weight(i), quadrature->weight(i))) << "Unexpected weight in copied quadrature.";
                expect(eq(copied_quadrature->node(i)[0], quadrature->node(i)[0])) << "Unexpected node coordinate 0 in copied quadrature.";
                expect(eq(copied_quadrature->node(i)[1], quadrature->node(i)[1])) << "Unexpected node coordinate 1 in copied quadrature.";
            }
        };
    }
};

}
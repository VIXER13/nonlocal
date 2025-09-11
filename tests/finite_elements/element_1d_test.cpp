#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

template<size_t Element_Order, size_t Quadrature_Order>
std::unique_ptr<element_1d_integrate_base<T>> make_element() {
    return std::make_unique<element_1d_integrate<T, lagrangian_element_1d, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

std::vector<std::unique_ptr<element_1d_integrate_base<T>>> init_lagrangian_elements_1d() {
    std::vector<std::unique_ptr<element_1d_integrate_base<T>>> result;
    result.emplace_back(make_element<0, 1>());
    result.emplace_back(make_element<1, 1>());
    result.emplace_back(make_element<2, 2>());
    result.emplace_back(make_element<3, 2>());
    result.emplace_back(make_element<4, 3>());
    result.emplace_back(make_element<5, 3>());
    return result;
}

struct basis_summator final {
    const element_1d_integrate_base<T>& element;
    T point = T{0};

    std::pair<T, T> operator()(const std::pair<T, T>& sum, const size_t i) const {
        return {
            sum.first + element.N(i, point),
            sum.second + element.Nxi(i, point)
        };
    }
};

const suite<"element_1d"> _ = [] {
    static constexpr std::array<size_t, 6> Qnodes_Count = {1, 1, 2, 2, 3, 3};
    size_t order = 0;
    for(const auto& element : init_lagrangian_elements_1d()) {
        const std::string suffix = "_order_" + std::to_string(order);

        test("nodes_count" + suffix) = [&element, order] {
            expect(eq(element->nodes_count(), order + 1)) << "Unexpected nodes count.";
        };

        test("qnodes_count" + suffix) = [&element, count = Qnodes_Count[order]] {
            expect(eq(element->qnodes_count(), count)) << "Unexpected qnodes count.";
        };

        test("boundaries" + suffix) = [&element] {
            expect(eq(element->boundary(side_1d::LEFT), T{-1})) << "Unexpected left boundary position.";
            expect(eq(element->boundary(side_1d::RIGHT), T{1})) << "Unexpected right boundary position.";
        };

        test("node_values" + suffix) = [&element] {
            for(const size_t i : element->nodes())
                for(const size_t j : element->nodes())
                    expect(eq(element->N(i, element->node(j)), T(i == j))) <<
                        "Unexpected value of function " + std::to_string(i) + " at node " + std::to_string(j);
        };

        const T element_length = element->boundary(side_1d::RIGHT) - element->boundary(side_1d::LEFT);

        test("weights_sum" + suffix) = [&element, element_length, qnodes = element->qnodes()] {
            const auto weight_summator = [&element](const T sum, const size_t qnode) { return sum + element->weight(qnode); };
            const T weights_sum = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, weight_summator);
            expect(lt(std::abs(weights_sum - element_length), std::numeric_limits<T>::epsilon())) <<
                "The weights sum does not match the element length.";
        };

        test("element_integral" + suffix) = [&element, element_length, qnodes = element->qnodes()] {
            const auto integrator = [&element](const T sum, const size_t q) {
                const auto nodes = element->nodes();
                const auto summator = [&element, q](const T sum, const size_t i) { return sum + element->qN(i, q); };
                return sum + element->weight(q) * std::accumulate(nodes.begin(), nodes.end(), T{0}, summator);
            };
            const T integral = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, integrator);
            static constexpr T Epsilon = T{5e-16};
            expect(lt(std::abs(integral - element_length), Epsilon)) << 
                "The sum of the integrals of all basis functions does not match with the element length.";
        };
        
        test("basis_sum" + suffix) = [&element, nodes = element->nodes()] {
            for(const T point : {T{-0.89}, T{-0.5234}, T{-0.03}, T{0.11}, T{0.57}, T{0.9844}}) { // some test points
                const basis_summator summator{*element, point};
                const auto [N_sum, Nxi_sum] = std::accumulate(nodes.begin(), nodes.end(), std::pair{T{0}, T{0}}, summator);
                static constexpr T Epsilon_Basis = T{4e-16};
                expect(lt(std::abs(N_sum - T{1}), Epsilon_Basis)) << 
                    "Unexpected basis functions sum in point = " + std::to_string(point);
                static constexpr T Epsilon_Derivative = T{5e-15};
                expect(lt(std::abs(Nxi_sum), Epsilon_Derivative)) << 
                    "Unexpected basis functions derivatives sum in point = " + std::to_string(point);
            }                
        };

        ++order;
    }
};

}
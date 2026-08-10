#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <string>
#include <ranges>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

constexpr auto Sample_Nodes = std::array{T{-0.89}, T{-0.5234}, T{-0.03}, T{0.11}, T{0.57}, T{0.9844}};

template<size_t Order>
std::unique_ptr<element_1d_base<T>> make_element_1d() {
    return std::make_unique<element_1d<T, lagrangian_element_1d, Order>>();
}

template<size_t Order>
std::unique_ptr<quadrature_1d_base<T>> make_quadrature_1d() {
    return std::make_unique<quadrature_1d<T, gauss, Order>>();
}

std::unique_ptr<element_1d_base<T>> make_element_1d(const size_t order) {
    switch(order) {
        case 1: return make_element_1d<1>();
        case 2: return make_element_1d<2>();
        case 3: return make_element_1d<3>();
        case 4: return make_element_1d<4>();
        case 5: return make_element_1d<5>();
        default: throw std::domain_error{"Unsupported element order " + std::to_string(order)};
    }
}

std::unique_ptr<quadrature_1d_base<T>> make_quadrature_1d(const size_t order) {
    switch(order) {
        case 1: return make_quadrature_1d<1>();
        case 2: return make_quadrature_1d<2>();
        case 3: return make_quadrature_1d<3>();
        case 4: return make_quadrature_1d<4>();
        case 5: return make_quadrature_1d<5>();
        default: throw std::domain_error{"Unsupported quadrature order " + std::to_string(order)};
    }
}

const suite<"element_1d"> _ = [] {
    for(const size_t element_order : std::ranges::iota_view(1zu, 6zu)) {
        const std::string suffix = "_element_order_" + std::to_string(element_order);
        const auto element = make_element_1d(element_order);

        test("nodes_count" + suffix) = [&element, element_order] {
            expect(eq(element->nodes_count(), element_order + 1)) << "Unexpected nodes count.";
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

        test("basis_sum" + suffix) = [&element] {
            const auto nodes = element->nodes();
            for(const T point : Sample_Nodes) {
                const auto basis_sum = std::accumulate(nodes.begin(), nodes.end(), T{0},
                    [&element, point](const T sum, const size_t i) { return sum + element->N(i, point); });
                static constexpr T Epsilon_Basis = T{4e-16};
                expect(approx(basis_sum, T{1}, Epsilon_Basis)) << 
                    "Unexpected basis functions sum in point = " + std::to_string(point);
            }                
        };

        test("basis_derivative_sum" + suffix) = [&element] {
            const auto nodes = element->nodes();
            for(const T point : Sample_Nodes) {
                const auto basis_derivative_sum = std::accumulate(nodes.begin(), nodes.end(), T{0},
                    [&element, point](const T sum, const size_t i) { return sum + element->Nxi(i, point); });
                static constexpr T Epsilon_Derivative = T{5e-15};
                expect(approx(basis_derivative_sum, T{0}, Epsilon_Derivative)) << 
                    "Unexpected basis functions derivatives sum in point = " + std::to_string(point);
            }                
        };

        test("copy" + suffix) = [&element] {
            const auto copy = element->copy();
            expect(neq(copy.get(), nullptr)) << "Copied element is nullptr.";
            expect(eq(copy->nodes_count(), element->nodes_count())) << "Unexpected nodes count in copied element.";
            expect(eq(copy->boundary(side_1d::LEFT), element->boundary(side_1d::LEFT))) << 
                "Unexpected left boundary position in copied element.";
            expect(eq(copy->boundary(side_1d::RIGHT), element->boundary(side_1d::RIGHT))) << 
                "Unexpected right boundary position in copied element.";
            for(const T node : Sample_Nodes) {
                for(const size_t i : element->nodes())
                    expect(eq(copy->N(i, node), element->N(i, node))) << 
                        "Unexpected value of function " + std::to_string(i) + " at node " + std::to_string(node) + " in copied element.";
            }
        };

        for (const size_t quadrature_order : std::ranges::iota_view{element_order, 6zu}) {
            const std::string quadrature_suffix = suffix + "_quadrature_order_" + std::to_string(quadrature_order);
            const element_1d_integrate<T> integrated_element{*element, *make_quadrature_1d(quadrature_order)};

            test("weights" + quadrature_suffix) = [&integrated_element] {
                const T element_length = integrated_element.element().boundary(side_1d::RIGHT) - integrated_element.element().boundary(side_1d::LEFT);
                const auto qnodes = integrated_element.qnodes();
                const T weights_sum = std::accumulate(qnodes.begin(), qnodes.end(), T{0},
                    [&integrated_element](const T sum, const size_t qnode) { return sum + integrated_element.weight(qnode); });
                expect(approx(weights_sum, element_length, std::numeric_limits<T>::epsilon())) << 
                    "The weights sum does not match the element length.";
            };

            test("integral" + quadrature_suffix) = [&integrated_element] {
                const T element_length = integrated_element.element().boundary(side_1d::RIGHT) - integrated_element.element().boundary(side_1d::LEFT);
                const auto qnodes = integrated_element.qnodes();
                const auto integrator = [&integrated_element](const T sum, const size_t q) {
                    const auto nodes = integrated_element.nodes();
                    const auto summator = [&integrated_element, q](const T sum, const size_t i) { return sum + integrated_element.qN(i, q); };
                    return sum + integrated_element.weight(q) * std::accumulate(nodes.begin(), nodes.end(), T{0}, summator);
                };
                const T integral = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, integrator);
                static constexpr T Epsilon = T{5e-16};
                expect(approx(integral, element_length, Epsilon)) << 
                    "The sum of the integrals of all basis functions does not match with the element length.";
            };

            test("copy" + quadrature_suffix) = [&integrated_element] {
                const auto copy_ptr = integrated_element.copy();
                const auto& copy = dynamic_cast<const element_1d_integrate<T>&>(*copy_ptr);
                expect(neq(&copy, nullptr)) << "Copied element is not of type element_1d_integrate.";
                expect(eq(copy.nodes_count(), integrated_element.nodes_count())) << "Unexpected nodes count in copied element.";
                expect(eq(copy.qnodes_count(), integrated_element.qnodes_count())) << "Unexpected qnodes count in copied element.";
                for(const size_t i : integrated_element.nodes())
                    expect(eq(copy.nearest_qnode(i), integrated_element.nearest_qnode(i))) << 
                        "Unexpected nearest qnode for node " + std::to_string(i) + " in copied element.";
                for(const size_t q : integrated_element.qnodes())
                    expect(eq(copy.weight(q), integrated_element.weight(q))) << 
                        "Unexpected weight for qnode " + std::to_string(q) + " in copied element.";
                for(const size_t i : integrated_element.nodes())
                    for(const size_t q : integrated_element.qnodes()) {
                        expect(eq(copy.qN(i, q), integrated_element.qN(i, q))) << 
                            "Unexpected value of function " + std::to_string(i) + " at qnode " + std::to_string(q) + " in copied element.";
                        expect(eq(copy.qNxi(i, q), integrated_element.qNxi(i, q))) << 
                            "Unexpected derivative of function " + std::to_string(i) + " at qnode " + std::to_string(q) + " in copied element.";
                    };
            };
        }
    }
};

}
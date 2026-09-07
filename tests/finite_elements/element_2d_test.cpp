#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>
#include <cstddef>
#include <memory>
#include <ranges>
#include <type_traits>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

constexpr std::array<std::array<T, 2>, 5> Sample_Rectangle_Nodes = {{
    {T{-0.89}, T{-0.57}},
    {T{0.11}, T{0.9844}},
    {T{-0.5234}, T{0.11}},
    {T{-0.03}, T{-0.57}},
    {T{0.57}, T{0.9844}}
}};

constexpr std::array<std::array<T, 2>, 5> Sample_Triangle_Nodes = {{
    {T{0.146}, T{0.34}},
    {T{0.732}, T{0.11}},
    {T{0.45}, T{0.37}},
    {T{0.25}, T{0.25}},
    {T{0.5}, T{0.1}}
}};

const suite<"element_2d"> _ = [] {
    for(const element_t type : {element_t::Triangle, element_t::Serendipity, element_t::Lagrangian}) {
        for(const size_t order : std::ranges::iota_view{1zu, 6zu}) {
            if (type == element_t::Triangle && order > 3)
                continue;
            const std::string type_name = type == element_t::Triangle ? "_triangle" :
                                          type == element_t::Serendipity ? "_serendipity" : "_lagrangian";
            const std::string suffix = type_name + "_order_" + std::to_string(order);
            const auto& sample_nodes = type == element_t::Triangle ? Sample_Triangle_Nodes : Sample_Rectangle_Nodes;
            const auto element = make_element_2d<T>(type, order);

            test("node_values" + suffix) = [&element] {
                static constexpr T Epsilon = T{1e-15};
                for(const size_t i : element->nodes())
                    for(const size_t j : element->nodes())
                        expect(approx(element->N(i, element->node(j)), T(i == j), Epsilon)) <<
                            "Unexpected value of function " + std::to_string(i) + " at node " + std::to_string(j);
            };

            test("basis_sum" + suffix) = [&element, &sample_nodes] {
                const auto nodes = element->nodes();
                for(const auto& point : sample_nodes) {
                    const auto basis_sum = std::accumulate(nodes.begin(), nodes.end(), T{0},
                        [&element, point](const T sum, const size_t i) { return sum + element->N(i, point); });
                    static constexpr T Epsilon_Basis = T{2e-15};
                    expect(approx(basis_sum, T{1}, Epsilon_Basis)) << 
                        "Unexpected basis functions sum in point = (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ")";
                }                
            };

            test("basis_derivative_xi_sum" + suffix) = [&element, &sample_nodes] {
                const auto nodes = element->nodes();
                for(const auto& point : sample_nodes) {
                    const auto basis_derivative_sum = std::accumulate(nodes.begin(), nodes.end(), T{0},
                        [&element, point](const T sum, const size_t i) { return sum + element->Nxi(i, point); });
                    static constexpr T Epsilon_Derivative = T{4e-15};
                    expect(approx(basis_derivative_sum, T{0}, Epsilon_Derivative)) << 
                        "Unexpected basis functions derivatives sum in point = (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ")";
                }                
            };

            test("basis_derivative_eta_sum" + suffix) = [&element, &sample_nodes] {
                const auto nodes = element->nodes();
                for(const auto& point : sample_nodes) {
                    const auto basis_derivative_sum = std::accumulate(nodes.begin(), nodes.end(), T{0},
                        [&element, point](const T sum, const size_t i) { return sum + element->Neta(i, point); });
                    static constexpr T Epsilon_Derivative = T{4e-15};
                    expect(approx(basis_derivative_sum, T{0}, Epsilon_Derivative)) << 
                        "Unexpected basis functions derivatives sum in point = (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ")";
                }                
            };

            test("copy" + suffix) = [&element, &sample_nodes] {
                const auto copy = element->copy();
                expect(neq(copy.get(), nullptr)) << "Copied element is nullptr.";
                expect(eq(copy->nodes_count(), element->nodes_count())) << "Unexpected nodes count in copied element.";
                for(const size_t i : element->nodes()) {
                    expect(eq(copy->node(i)[0], element->node(i)[0])) << "Unexpected node position in copied element.";
                    expect(eq(copy->node(i)[1], element->node(i)[1])) << "Unexpected node position in copied element.";
                    for(const auto& point : sample_nodes) {
                        expect(eq(copy->N(i, point), element->N(i, point))) << 
                            "Unexpected value of function " + std::to_string(i) + " at point (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ") in copied element.";
                        expect(eq(copy->Nxi(i, point), element->Nxi(i, point))) << 
                            "Unexpected value of function derivative " + std::to_string(i) + " at point (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ") in copied element.";
                        expect(eq(copy->Neta(i, point), element->Neta(i, point))) << 
                            "Unexpected value of function derivative " + std::to_string(i) + " at point (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ") in copied element.";
                    }
                }
            };

            for(const size_t quadrature_order : std::ranges::iota_view{order, 6zu}) {
                const std::string quadrature_suffix = suffix + "_quadrature_order_" + std::to_string(quadrature_order);
                const auto geometry_type = type == element_t::Triangle ? geometry_t::Triangle : geometry_t::Rectangle;
                const element_2d_integrate<T> integrated_element{*element, *make_quadrature_2d<T>(geometry_type, quadrature_order)};

                test("weights_sum" + quadrature_suffix) = [&integrated_element] {
                    const auto qnodes = integrated_element.qnodes();
                    const auto weight_summator = [&integrated_element](const T sum, const size_t qnode) {
                        return sum + integrated_element.weight(qnode);
                    };
                    const T weights_sum = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, weight_summator);
                    static constexpr T Epsilon = 1e-15;
                    expect(approx(weights_sum, dynamic_cast<const geometry_2d_base<T>&>(integrated_element.element()).area(), Epsilon)) <<
                        "Unexpected weights sum. The sum of the integrals of all basis functions does not match with the element area.";
                };

                test("integral" + quadrature_suffix) = [&integrated_element] {
                    const auto qnodes = integrated_element.qnodes();
                    const auto integrator = [&integrated_element](const T sum, const size_t q) {
                        const auto nodes = integrated_element.nodes();
                        const auto summator = [&integrated_element, q](const T sum, const size_t i) { return sum + integrated_element.qN(i, q); };
                        return sum + integrated_element.weight(q) * std::accumulate(nodes.begin(), nodes.end(), T{0}, summator);
                    };
                    const T integral = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, integrator);
                    static constexpr T Epsilon = 5e-15;
                    expect(approx(integral, dynamic_cast<const geometry_2d_base<T>&>(integrated_element.element()).area(), Epsilon)) <<
                        "The sum of the integrals of all basis functions does not match with the element area.";
                };

                test("integrated_copy" + quadrature_suffix) = [&integrated_element] {
                    const auto copy_ptr = integrated_element.copy();
                    const auto& copy = dynamic_cast<const element_2d_integrate<T>&>(*copy_ptr);
                    expect(neq(&copy, nullptr)) << "Copied element is not of type element_2d_integrate.";
                    expect(eq(copy.nodes_count(), integrated_element.nodes_count())) << "Unexpected nodes count in copied element.";
                    expect(eq(copy.qnodes_count(), integrated_element.qnodes_count())) << "Unexpected quadrature nodes count in copied element.";
                    for(const size_t i : integrated_element.nodes()) {
                        expect(eq(copy.nearest_qnode(i), integrated_element.nearest_qnode(i))) << "Unexpected nearest quadrature node in copied element.";
                        expect(eq(copy.element().node(i)[0], integrated_element.element().node(i)[0])) << "Unexpected node position in copied element.";
                        expect(eq(copy.element().node(i)[1], integrated_element.element().node(i)[1])) << "Unexpected node position in copied element.";
                    }
                    for(const size_t q : integrated_element.qnodes())
                        expect(eq(copy.weight(q), integrated_element.weight(q))) << "Unexpected quadrature weight in copied element.";
                    for(const size_t i : integrated_element.nodes())
                        for(const size_t q : integrated_element.qnodes()) {
                            expect(eq(copy.qN(i, q), integrated_element.qN(i, q))) << "Unexpected value of function " + std::to_string(i) + " at quadrature node " + std::to_string(q) + " in copied element.";
                            expect(eq(copy.qNxi(i, q), integrated_element.qNxi(i, q))) << "Unexpected value of function derivative " + std::to_string(i) + " at quadrature node " + std::to_string(q) + " in copied element.";
                            expect(eq(copy.qNeta(i, q), integrated_element.qNeta(i, q))) << "Unexpected value of function derivative " + std::to_string(i) + " at quadrature node " + std::to_string(q) + " in copied element.";
                        }
                };
            }
        }
    }
};

}
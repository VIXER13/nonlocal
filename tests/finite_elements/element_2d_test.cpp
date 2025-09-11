#include <metamath/metamath.hpp>

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;
using T = double;

template<template<class, auto...> class Element_Type, size_t Quadrature_Order_X1, size_t Quadrature_Order_X2, size_t... Element_Order>
std::unique_ptr<element_2d_integrate_base<T>> make_element() {
    return std::make_unique<element_2d_integrate<T, Element_Type, Element_Order...>>(
        quadrature_1d<T, gauss, Quadrature_Order_X1>{},
        quadrature_1d<T, gauss, Quadrature_Order_X2>{}
    );
}

std::vector<std::unique_ptr<element_2d_integrate_base<T>>> init_lagrangian_elements_2d() {
    std::vector<std::unique_ptr<element_2d_integrate_base<T>>> result;
    result.emplace_back(make_element<lagrangian_element_2d, 1, 1, 0, 0>());
    result.emplace_back(make_element<lagrangian_element_2d, 1, 1, 0, 1>());
    result.emplace_back(make_element<lagrangian_element_2d, 1, 1, 1, 0>());
    result.emplace_back(make_element<lagrangian_element_2d, 1, 1, 1, 1>());
    result.emplace_back(make_element<lagrangian_element_2d, 1, 2, 1, 2>());
    result.emplace_back(make_element<lagrangian_element_2d, 2, 1, 2, 1>());
    result.emplace_back(make_element<lagrangian_element_2d, 2, 2, 2, 2>());
    result.emplace_back(make_element<lagrangian_element_2d, 2, 2, 2, 3>());
    result.emplace_back(make_element<lagrangian_element_2d, 2, 2, 3, 2>());
    result.emplace_back(make_element<lagrangian_element_2d, 2, 2, 3, 3>());
    return result;
}

std::vector<std::unique_ptr<element_2d_integrate_base<T>>> init_serendipity_elements_2d() {
    std::vector<std::unique_ptr<element_2d_integrate_base<T>>> result;
    result.emplace_back(make_element<serendipity, 1, 1, 0>());
    result.emplace_back(make_element<serendipity, 1, 1, 1>());
    // Second and third order shall be created manually, due to a bug in compiling specialized templates
    result.emplace_back(std::make_unique<element_2d_integrate<T, serendipity, 2>>(quadrature_1d<T, gauss, 2>{}, quadrature_1d<T, gauss, 2>{}));
    result.emplace_back(std::make_unique<element_2d_integrate<T, serendipity, 3>>(quadrature_1d<T, gauss, 2>{}, quadrature_1d<T, gauss, 2>{}));
    result.emplace_back(make_element<serendipity, 3, 3, 4>());
    result.emplace_back(make_element<serendipity, 3, 3, 5>());
    return result;
}

std::vector<std::unique_ptr<element_2d_integrate_base<T>>> init_triangle_elements_2d() {
    std::vector<std::unique_ptr<element_2d_integrate_base<T>>> result;
    result.emplace_back(make_element<triangle, 1, 1, 0>());
    result.emplace_back(make_element<triangle, 1, 1, 1>());
    result.emplace_back(make_element<triangle, 2, 2, 2>());
    result.emplace_back(make_element<triangle, 2, 2, 3>());
    return result;
}

struct basis_summator final {
    const element_2d_integrate_base<T>& element;
    std::array<T, 2> point = {};

    std::array<T, 3> operator()(const std::array<T, 3>& sum, const size_t i) const {
        return {
            sum[0] + element.N(i, point),
            sum[1] + element.Nxi(i, point),
            sum[2] + element.Neta(i, point)
        };
    }
};

std::string to_string(const std::array<T, 2>& point) {
    return '(' + std::to_string(point[0]) + " ," + std::to_string(point[1]) + ')';
}

std::vector<std::array<T, 2>> init_points(const std::string& element_type) {
    std::vector<std::array<T, 2>> points;
    if (element_type == "triangle") {
        points = {
            {T{0.146}, T{0.34}},
            {T{0.732}, T{0.11}},
            {T{0.45}, T{0.37}}
        };
    } else {
        points.reserve(9);
        for(const T x : {T{-0.89}, T{0.5234}, T{-0.03}}) // some test points
            for(const T y : {T{0.11}, T{-0.57}, T{0.9844}})
                points.push_back({x, y});
    }
    return points;
}

const suite<"element_2d"> _ = [] {
    using namespace std::literals;
    for(const auto element_type : {"triangle"s, "serendipity"s, "lagrangian"s}) {
        size_t order = 0;
        const T element_area = element_type == "triangle" ? T{0.5} : T{4};
        const std::vector<std::array<T, 2>> points = init_points(element_type);
        for(const auto& element : element_type == "triangle"    ? init_triangle_elements_2d() :
                                  element_type == "serendipity" ? init_serendipity_elements_2d() :
                                                                  init_lagrangian_elements_2d()) {
            const std::string suffix = "_" + element_type + "_order_" + std::to_string(order);

            test("node_values" + suffix) = [&element] {
                static constexpr T Epsilon = std::is_same_v<T, float> ? T{1e-6} : T{1e-15};
                for(const size_t i : element->nodes())
                    for(const size_t j : element->nodes())
                        expect(lt(std::abs(element->N(i, element->node(j)) - T(i == j)), Epsilon)) <<
                            "Unexpected value of function " + std::to_string(i) + " at node " + std::to_string(j);
            };

            test("weights_sum" + suffix) = [&element, element_area, qnodes = element->qnodes()] {
                const auto weight_summator = [&element](const T sum, const size_t qnode) { return sum + element->weight(qnode); };
                const T weights_sum = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, weight_summator);
                static constexpr T Epsilon = T{1e-16};
                expect(lt(std::abs(weights_sum - element_area), Epsilon)) <<
                    "The weights sum does not match the element area.";
            };

            test("element_integral" + suffix) = [&element, element_area] {
                const auto integrator = [&element](const T sum, const size_t q) {
                    const auto nodes = element->nodes();
                    const auto& summator = [&element, q](const T sum, const size_t i) { return sum + element->qN(i, q); };
                    return sum + element->weight(q) * std::accumulate(nodes.begin(), nodes.end(), T{0}, summator);
                };
                const auto qnodes = element->qnodes();
                const T integral = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, integrator);
                static constexpr T Epsilon = T{1e-15};
                expect(lt(std::abs(integral - element_area), Epsilon)) << 
                    "The sum of the integrals of all basis functions does not match with the element area.";
            };

            basis_summator summator{*element};
            test("basis_sum" + suffix) = [&element, &summator, &points, nodes = element->nodes()] {
                for(const std::array<T, 2>& point : points) {
                    summator.point = point;
                    const auto [N_sum, Nxi_sum, Neta_sum] = std::accumulate(nodes.begin(), nodes.end(), std::array{T{0}, T{0}, T{0}}, summator);
                    static constexpr T Epsilon_Basis = T{1e-15};
                    expect(lt(std::abs(N_sum - T{1}), Epsilon_Basis)) << 
                        "Unexpected basis functions sum in point = " + to_string(point);
                    static constexpr T Epsilon_Derivative = T{3e-15};
                    expect(lt(std::abs(Nxi_sum), Epsilon_Derivative)) << 
                        "Unexpected Nxi sum in point = " + to_string(point);
                    expect(lt(std::abs(Neta_sum), Epsilon_Derivative)) << 
                        "Unexpected Neta sum in point = " + to_string(point);
                }
            };

            ++order;
        }
    }
};

}
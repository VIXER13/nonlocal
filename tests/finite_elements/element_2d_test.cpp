#include "init_elements.hpp"

#include <boost/ut.hpp>

#include <numeric>

namespace {

using namespace boost::ut;
using namespace metamath::finite_element;

template<class T>
struct basis_summator final {
    const element_2d_integrate_base<T>& element;
    std::array<T, 2> point = {};

    std::array<T, 3> operator()(const std::array<T, 3>& sum, const size_t i) {
        return {
            sum[0] + element.N(i, point),
            sum[1] + element.Nxi(i, point),
            sum[2] + element.Neta(i, point)
        };
    }
};

template<class T>
std::string to_string(const std::array<T, 2>& point) {
    return '(' + std::to_string(point[0]) + " ," + std::to_string(point[1]) + ')';
}

template<class T>
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

const suite _ = [] {
    using namespace std::literals;
    for(const auto element_type : {"triangle"s, "serendipity"s, "lagrangian"s}) {
        test(element_type + "_element_2d") = [element_type]<class T> {
            size_t order = 0;
            const T element_area = element_type == "triangle" ? T{0.5} : T{4};
            const std::vector<std::array<T, 2>> points = init_points<T>(element_type);

            for(const auto& element : element_type == "triangle" ? unit_tests::init_triangle_elements_2d<T>() :
                                      element_type == "serendipity" ? unit_tests::init_serendipity_elements_2d<T>() :
                                      unit_tests::init_lagrangian_elements_2d<T>()) {
                const std::string suffix = "_order_" + std::to_string(order) + '_' + std::string{reflection::type_name<T>()};

                test("node_values" + suffix) = [&element, nodes = element->nodes()] {
                    static constexpr T epsilon = std::is_same_v<T, float> ? T{1e-6} : T{1e-15};
                    for(const size_t i : nodes)
                        for(const size_t j : nodes)
                            expect(lt(std::abs(element->N(i, element->node(j)) - T(i == j)), epsilon)) <<
                                "Unexpected value of function " + std::to_string(i) + " at node " + std::to_string(j) + '.';
                };

                test("weights_sum" + suffix) = [&element, element_area] {
                    const auto qnodes = element->qnodes();
                    const auto weight_summator = [&element](const T sum, const size_t qnode) { return sum + element->weight(qnode); };
                    const T weights_sum = std::accumulate(qnodes.begin(), qnodes.end(), T{0}, weight_summator);
                    static constexpr T epsilon = std::is_same_v<T, float> ? T{1e-7} : T{1e-16};
                    expect(lt(std::abs(weights_sum - element_area), epsilon)) <<
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
                    static constexpr T epsilon = std::is_same_v<T, float> ? T{1e-6} : T{1e-15};
                    expect(lt(std::abs(integral - element_area), epsilon)) << 
                        "The sum of the integrals of all basis functions does not match with the element area.";
                };

                basis_summator<T> summator{*element};
                test("basis_sum" + suffix) = [&element, &summator, &points, nodes = element->nodes()] {
                    for(const std::array<T, 2>& point : points) {
                        summator.point = point;
                        const auto [N_sum, Nxi_sum, Neta_sum] = std::accumulate(nodes.begin(), nodes.end(), std::array{T{0}, T{0}, T{0}}, summator);
                        static constexpr T epsilon_N = std::is_same_v<T, float> ? T{1e-6} : T{1e-15};
                        expect(lt(std::abs(N_sum - T{1}), epsilon_N)) << 
                            "Unexpected basis functions sum in point = " + to_string(point) + '.';
                        static constexpr T epsilon_Nxi = std::is_same_v<T, float> ? T{2e-6} : T{3e-15};
                        expect(lt(std::abs(Nxi_sum), epsilon_Nxi)) << 
                            "Unexpected Nxi sum in point = " + to_string(point) + '.';
                        static constexpr T epsilon_Neta = std::is_same_v<T, float> ? T{2e-6} : T{3e-15};
                        expect(lt(std::abs(Neta_sum), epsilon_Neta)) << 
                            "Unexpected Neta sum in point = " + to_string(point) + '.';
                    }
                };

                test("basis_incomplete_sum" + suffix) = [&element, &summator, &points, nodes = element->nodes()] {
                    for(const std::array<T, 2>& point : points) {
                        summator.point = point;
                        const auto [N_sum, Nxi_sum, Neta_sum] = std::accumulate(nodes.begin(), std::ranges::prev(nodes.end()), std::array{T{0}, T{0}, T{0}}, summator);
                        static constexpr T epsilon_basis = T{1e-4};
                        expect(ge(std::abs(N_sum - T{1}), element->nodes_count() > 1 ? epsilon_basis : T{0})) << 
                            "The partial sum at point = " + to_string(point) + " must not be equal to 1.";
                        static constexpr T epsilon_derivatives = T{1e-3};
                        expect(ge(std::abs(Nxi_sum), element->nodes_count() > 2 ? epsilon_derivatives : T{0})) << 
                            "The partial sum of Nxi at point = " + to_string(point) + " must not be equal to 1.";
                        expect(ge(std::abs(Neta_sum), element->nodes_count() > 2 ? epsilon_derivatives : T{0})) << 
                            "The partial sum of Neta at point = " + to_string(point) + " must not be equal to 1.";
                    }
                };

                ++order;
            }
        } | std::tuple<double>{};
    }
};

}
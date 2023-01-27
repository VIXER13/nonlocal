#include "metamath.hpp"

#include <gtest/gtest.h>

#include <array>
#include <numeric>
#include <string_view>

namespace {

using namespace metamath::finite_element;

template<class T>
struct test_data_1d final {
    std::shared_ptr<element_1d_integrate_base<T>> element;
    size_t expected_nodes_count;
    size_t expected_qnodes_count;
    T epsilon = std::numeric_limits<T>::epsilon();
};

template<class T>
struct test_data_2d final {
    std::shared_ptr<element_2d_integrate_base<T>> element;
    size_t expected_nodes_count;
    size_t expected_qnodes_count;
    T expected_area;
    T epsilon = std::numeric_limits<T>::epsilon();
};

class elements_1d : public ::testing::TestWithParam<test_data_1d<float>> {};

class elements_2d : public ::testing::TestWithParam<test_data_2d<float>> {};

constexpr auto tests_names_2d = [](const ::testing::TestParamInfo<elements_2d::ParamType>& info) {
    using namespace std::literals;
    static constexpr auto names = std::array{
        "triangle_0"sv,
        "triangle_1"sv,
        "triangle_2"sv,
        "triangle_3"sv,
        "serendipity_0"sv,
        "serendipity_1"sv,
        "serendipity_2"sv,
        "serendipity_3"sv,
        "serendipity_4"sv,
        "serendipity_5"sv,
        "lagrange_0_0"sv,
        "lagrange_0_1"sv,
        "lagrange_1_0"sv,
        "lagrange_1_1"sv,
        "lagrange_1_2"sv,
        "lagrange_2_1"sv,
        "lagrange_2_2"sv,
        "lagrange_2_3"sv,
        "lagrange_3_2"sv,
        "lagrange_3_3"sv
    };
    return names[info.index].data();
};

template<class T, size_t Element_Order, size_t Quadrature_Order>
std::shared_ptr<element_1d_integrate_base<T>> make_element_1d() {
    return std::make_shared<element_1d_integrate<T, lagrangian_element_1d, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

template<class T, template<class, size_t> class Element_Kind, size_t Element_Order, size_t Quadrature_Order>
std::shared_ptr<element_2d_integrate_base<T>> make_element_2d() {
    return std::make_shared<element_2d_integrate<T, Element_Kind, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{},
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

template<class T, size_t Element_Order_X, size_t Element_Order_Y, size_t Quadrature_Order_X, size_t Quadrature_Order_Y>
std::shared_ptr<element_2d_integrate_base<T>> make_lagrangian_element_2d() {
    return std::make_shared<element_2d_integrate<T, lagrangian_element_2d, Element_Order_X, Element_Order_Y>>(
        quadrature_1d<T, gauss, Quadrature_Order_X>{},
        quadrature_1d<T, gauss, Quadrature_Order_Y>{}
    );
}

template<class T>
T weights_sum(const element_integrate_base<T>& element) {
    const auto weight_summator = [&element](const T sum, const size_t q) { return sum + element.weight(q); };
    return std::accumulate(element.qnodes().begin(), element.qnodes().end(), T{0}, weight_summator);
}

template<class T>
T integrate_element(const element_integrate_base<T>& element) {
    const auto integrator = [&element](const T sum, const size_t q) {
        const auto& summator = [&element, q](const T sum, const size_t i) { return sum + element.qN(i, q); };
        return sum + element.weight(q) * std::accumulate(element.nodes().begin(), element.nodes().end(), T{0}, summator);
    };
    return std::accumulate(element.qnodes().begin(), element.qnodes().end(), T{0}, integrator);
}

}

TEST_P(elements_1d, lagrangian) {
    const auto& data = GetParam();

    ASSERT_NE(data.element.get(), nullptr);
    const auto& element = *data.element;
    EXPECT_EQ(element.nodes_count(), data.expected_nodes_count);
    EXPECT_EQ(element.qnodes_count(), data.expected_qnodes_count);

    for(const size_t i : element.nodes())
        for(const size_t j : element.nodes())
            EXPECT_NEAR(element.N(i, element.node(j)), float(i == j), data.epsilon);

    const float element_length = element.boundary(side_1d::RIGHT) - element.boundary(side_1d::LEFT);
    EXPECT_NEAR(weights_sum(element), element_length, data.epsilon);
    EXPECT_NEAR(integrate_element(element), element_length, data.epsilon);

    for(const float point : {-0.89f, -0.5234f, -0.03f, 0.11f, 0.57f, 0.9844f}) { // some test points
        const auto basis_summator = [&element, point](const std::pair<float, float> sum, const size_t i) {
            return std::pair{sum.first + element.N(i, point), sum.second + element.Nxi(i, point)};
        };

        const auto [N_sum, Nxi_sum] = std::accumulate(element.nodes().begin(), element.nodes().end(), std::pair{0.f, 0.f}, basis_summator);
        EXPECT_NEAR(N_sum, 1.f, data.epsilon);
        EXPECT_NEAR(Nxi_sum, 0.f, data.epsilon);

        const auto [N_incomplete_sum, Nxi_incomplete_sum] =
            std::accumulate(element.nodes().begin(), std::ranges::prev(element.nodes().end()), std::pair{0.f, 0.f}, basis_summator);
        EXPECT_GE(std::abs(N_incomplete_sum - 1.f), element.nodes_count() > 1 ? data.epsilon : 0.f);
        EXPECT_GE(std::abs(Nxi_incomplete_sum), element.nodes_count() > 1 ? data.epsilon : 0.f);
    }
}

INSTANTIATE_TEST_CASE_P(
    elements,
    elements_1d,
    ::testing::Values(
        test_data_1d{make_element_1d<float, 0, 1>(), 1, 1},
        test_data_1d{make_element_1d<float, 1, 1>(), 2, 1},
        test_data_1d{make_element_1d<float, 2, 2>(), 3, 2, 1e-6f},
        test_data_1d{make_element_1d<float, 3, 2>(), 4, 2, 1e-6f},
        test_data_1d{make_element_1d<float, 4, 3>(), 5, 3, 1e-6f},
        test_data_1d{make_element_1d<float, 5, 3>(), 6, 3, 2.5e-6f}
    )
);

TEST_P(elements_2d, basises_2d) {
    const auto& data = GetParam();

    ASSERT_NE(data.element.get(), nullptr);
    const auto& element = *data.element;
    EXPECT_EQ(element.nodes_count(), data.expected_nodes_count);
    EXPECT_EQ(element.qnodes_count(), data.expected_qnodes_count);

    for(const size_t i : element.nodes())
        for(const size_t j : element.nodes())
            EXPECT_NEAR(element.N(i, element.node(j)), float(i == j), data.epsilon);

    EXPECT_NEAR(weights_sum(element), data.expected_area, data.epsilon);
    EXPECT_NEAR(integrate_element(element), data.expected_area, data.epsilon);

    for(const float x : {-0.89f, -0.5234f, -0.03f}) // some test points
        for(const float y : {0.11f, 0.57f, 0.9844f}) {
            const auto basis_summator = [&element, point = std::array{x, y}](const std::array<float, 3>& sum, const size_t i) {
                return std::array{
                    sum[0] + element.N(i, point),
                    sum[1] + element.Nxi(i, point),
                    sum[2] + element.Neta(i, point)
                };
            };

            const auto [N_sum, Nxi_sum, Neta_sum] = std::accumulate(element.nodes().begin(), element.nodes().end(), std::array{0.f, 0.f, 0.f}, basis_summator);
            EXPECT_NEAR(N_sum, 1.f, data.epsilon);
            EXPECT_NEAR(Nxi_sum, 0.f, data.epsilon);
            EXPECT_NEAR(Neta_sum, 0.f, data.epsilon);

            const auto [N_incomplete_sum, Nxi_incomplete_sum, Neta_incomplete_sum] =
                std::accumulate(element.nodes().begin(), std::ranges::prev(element.nodes().end()), std::array{0.f, 0.f, 0.f}, basis_summator);
            EXPECT_GE(std::abs(N_incomplete_sum - 1.f), element.nodes_count() > 1 ? data.epsilon : 0.f);
            EXPECT_GE(std::abs(Nxi_incomplete_sum), element.nodes_count() > 2 ? data.epsilon : 0.f);
            EXPECT_GE(std::abs(Neta_incomplete_sum), element.nodes_count() > 2 ? data.epsilon : 0.f);
        }
}

INSTANTIATE_TEST_CASE_P(
    elements,
    elements_2d,
    ::testing::Values(
        test_data_2d{make_element_2d<float, triangle, 0, 1>(), 1, 1, 0.5f},
        test_data_2d{make_element_2d<float, triangle, 1, 1>(), 3, 1, 0.5f},
        test_data_2d{make_element_2d<float, triangle, 2, 2>(), 6, 4, 0.5f, 1e-6f},
        test_data_2d{make_element_2d<float, triangle, 3, 2>(), 10, 4, 0.5f, 1e-5f},
        test_data_2d{make_element_2d<float, serendipity, 0, 1>(), 1, 1, 4.f},
        test_data_2d{make_element_2d<float, serendipity, 1, 1>(), 4, 1, 4.f},
        // When I call the make_element_2d function, the compiler does not see the specialization for the 2nd and 3rd degree. Looks like a gcc 11.3 compiler bug.
        test_data_2d<float>{std::make_shared<element_2d_integrate<float, serendipity, 2>>(quadrature_1d<float, gauss, 2>{}, quadrature_1d<float, gauss, 2>{}), 8, 4, 4.f, 1e-6f},
        test_data_2d<float>{std::make_shared<element_2d_integrate<float, serendipity, 3>>(quadrature_1d<float, gauss, 2>{}, quadrature_1d<float, gauss, 2>{}), 12, 4, 4.f, 1e-6f},
        test_data_2d{make_element_2d<float, serendipity, 4, 3>(), 16, 9, 4.f, 2e-6f},
        test_data_2d{make_element_2d<float, serendipity, 5, 3>(), 20, 9, 4.f, 2e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 0, 0, 1, 1>(), 1, 1, 4.f},
        test_data_2d{make_lagrangian_element_2d<float, 0, 1, 1, 1>(), 2, 1, 4.f},
        test_data_2d{make_lagrangian_element_2d<float, 1, 0, 1, 1>(), 2, 1, 4.f},
        test_data_2d{make_lagrangian_element_2d<float, 1, 1, 1, 1>(), 4, 1, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 1, 2, 1, 2>(), 6, 2, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 2, 1, 2, 1>(), 6, 2, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 2, 2, 2, 2>(), 9, 4, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 2, 3, 2, 2>(), 12, 4, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 3, 2, 2, 2>(), 12, 4, 4.f, 1e-6f},
        test_data_2d{make_lagrangian_element_2d<float, 3, 3, 2, 2>(), 16, 4, 4.f, 1e-6f}
    ),
    tests_names_2d
);
#include "metamath.hpp"

#include <gtest/gtest.h>

#include <numeric>

namespace {

using namespace metamath::finite_element;

template<class T>
struct test_data final {
    std::shared_ptr<element_1d_base<T>> element;
    size_t expected_nodes_count;
    T epsilon;
};

class elements_1d : public ::testing::TestWithParam<test_data<float>> {};

}

TEST_P(elements_1d, lagrangian) {
    const auto& [element, expected_nodes_count, epsilon] = GetParam();

    ASSERT_NE(element.get(), nullptr);
    EXPECT_EQ(element->nodes_count(), expected_nodes_count);

    for(const size_t i : element->nodes())
        for(const size_t j : element->nodes())
            EXPECT_NEAR(element->N(i, element->node(j)), float(i == j), std::numeric_limits<float>::epsilon());

    for(const float point : {-0.89f, -0.5234f, -0.03f, 0.11f, 0.57f, 0.9844f}) { // some test points
        const auto basis_summator = [&element, point](const std::pair<float, float> sum, const size_t i) {
            return std::pair{sum.first + element->N(i, point), sum.second + element->Nxi(i, point)};
        };

        const auto [N_sum, Nxi_sum] = std::accumulate(element->nodes().begin(), element->nodes().end(), std::pair{0.f, 0.f}, basis_summator);
        EXPECT_NEAR(N_sum, 1.f, epsilon);
        EXPECT_NEAR(Nxi_sum, 0.f, epsilon);

        const auto [N_incomplete_sum, Nxi_incomplete_sum] =
            std::accumulate(element->nodes().begin(), std::prev(element->nodes().end()), std::pair{0.f, 0.f}, basis_summator);
        EXPECT_GE(std::abs(N_incomplete_sum - 1.f), element->nodes_count() > 1 ? epsilon : 0.f);
        EXPECT_GE(std::abs(Nxi_incomplete_sum), element->nodes_count() > 1 ? epsilon : 0.f);
    }
}

INSTANTIATE_TEST_CASE_P(
    basises,
    elements_1d,
    ::testing::Values(
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 0>>(), 1, std::numeric_limits<float>::epsilon()},
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 1>>(), 2, std::numeric_limits<float>::epsilon()},
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 2>>(), 3, 1e-6f},
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 3>>(), 4, 1e-6f},
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 4>>(), 5, 1e-6f},
        test_data<float>{std::make_shared<element_1d<float, lagrangian_element_1d, 5>>(), 6, 2.5e-6f}
    )
);
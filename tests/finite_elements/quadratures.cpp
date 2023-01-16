#include "metamath.hpp"

#include <gtest/gtest.h>

#include <numeric>

namespace {

using namespace metamath::finite_element;

template<class T>
struct test_data final {
    std::shared_ptr<quadrature_1d_base<T>> quadrature;
    size_t expected_nodes_count;
};

class gaussian_quadratures : public ::testing::TestWithParam<test_data<float>> {};

}

TEST_P(gaussian_quadratures, gauss) {
    const auto& [quadrature, expected_nodes_count] = GetParam();
    
    ASSERT_NE(quadrature.get(), nullptr);
    EXPECT_EQ(quadrature->nodes_count(), expected_nodes_count);

    EXPECT_EQ(quadrature->boundary(side_1d::LEFT), -1.f);
    EXPECT_EQ(quadrature->boundary(side_1d::RIGHT), 1.f);

    const auto weight_summator = [&quadrature](const float sum, const size_t node) { return sum + quadrature->weight(node); };
    const float weights_sum = std::accumulate(quadrature->nodes().begin(), quadrature->nodes().end(), 0.f, weight_summator);
    const float length = quadrature->boundary(side_1d::RIGHT) - quadrature->boundary(side_1d::LEFT);
    EXPECT_NEAR(weights_sum, length, std::numeric_limits<float>::epsilon());
}

INSTANTIATE_TEST_CASE_P(
    quadratures,
    gaussian_quadratures,
    ::testing::Values(
        test_data<float>{std::make_shared<quadrature_1d<float, gauss, 1>>(), 1},
        test_data<float>{std::make_shared<quadrature_1d<float, gauss, 2>>(), 2},
        test_data<float>{std::make_shared<quadrature_1d<float, gauss, 3>>(), 3},
        test_data<float>{std::make_shared<quadrature_1d<float, gauss, 4>>(), 4},
        test_data<float>{std::make_shared<quadrature_1d<float, gauss, 5>>(), 5} 
    ),
    [](const ::testing::TestParamInfo<gaussian_quadratures::ParamType>& info) {
        return std::to_string(info.index + 1);
    }
);
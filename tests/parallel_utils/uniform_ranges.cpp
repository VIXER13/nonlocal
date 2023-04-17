#include "init_uniform_ranges.hpp"

#include <gtest/gtest.h>

namespace {

using namespace parallel_utils;

constexpr size_t SIZE = 3;
constexpr size_t COUNTS = 6;

}

TEST(uniform_ranges, test_ranges_count) {
    EXPECT_ANY_THROW(init_uniform_ranges(SIZE, 0u));
    for(const size_t count : std::ranges::iota_view{1u, COUNTS}) {
        const auto ranges = init_uniform_ranges(SIZE, count);
        EXPECT_EQ(ranges.size(), count);
    }
}

TEST(uniform_ranges, one_range) {
    const auto one_range = init_uniform_ranges(SIZE, 1u);
    EXPECT_EQ(one_range[0].front(), 0u);
    EXPECT_EQ(one_range[0].back(),  SIZE - 1);
}

TEST(uniform_ranges, two_ranges) {
    const auto two_ranges = init_uniform_ranges(SIZE, 2u);
    EXPECT_EQ(two_ranges[0].front(), 0u);
    EXPECT_EQ(two_ranges[0].back(),  0u);
    EXPECT_EQ(two_ranges[1].front(), 1u);
    EXPECT_EQ(two_ranges[1].back(),  SIZE - 1);
}

TEST(uniform_ranges, three_ranges) {
    const auto three_ranges = init_uniform_ranges(SIZE, 3u);
    for(const size_t i : std::ranges::iota_view{0u, 3u}) {
        EXPECT_EQ(three_ranges[i].front(), i);
        EXPECT_EQ(three_ranges[i].back(),  i);
    }
}

TEST(uniform_ranges, four_ranges) {
    const auto four_ranges = init_uniform_ranges(SIZE, 4u);
    for(const size_t i : std::ranges::iota_view{0u, 3u}) {
        EXPECT_EQ(four_ranges[i].front(), i);
        EXPECT_EQ(four_ranges[i].back(),  i);
    }
    EXPECT_EQ(four_ranges[3].front(), SIZE);
    EXPECT_EQ(four_ranges[3].back(),  SIZE - 1);
}

TEST(uniform_ranges, more_ranges) {
    const auto more_ranges = init_uniform_ranges(SIZE, COUNTS);
    for(const size_t i : std::ranges::iota_view{0u, 3u}) {
        EXPECT_EQ(more_ranges[i].front(), i);
        EXPECT_EQ(more_ranges[i].back(),  i);
    }
    for(const size_t i : std::ranges::iota_view{3u, more_ranges.size()}) {
        EXPECT_EQ(more_ranges[i].front(), SIZE);
        EXPECT_EQ(more_ranges[i].back(),  SIZE - 1);
    }
}
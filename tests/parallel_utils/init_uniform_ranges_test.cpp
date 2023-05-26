#include "init_uniform_ranges.hpp"

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"uniform_ranges"> _ = [] {
    using namespace boost::ut;
    using namespace parallel_utils;

    static constexpr size_t SIZE = 3;
    static constexpr size_t COUNTS = 6;
    
    "ranges_count_0"_test = [] {
        expect(throws([] { init_uniform_ranges(SIZE, 0u); })) <<
            "When creating 0 ranges, an exception should be thrown.";
    };

    for(const size_t count : std::ranges::iota_view{1u, COUNTS})
        test("ranges_count_" + std::to_string(count)) = [count] {
            const auto ranges = init_uniform_ranges(SIZE, count);
            expect(eq(ranges.size(), count)) << "Unexpected ranges number.";
        };

    "range_1"_test = [] {
        const auto one_range = init_uniform_ranges(SIZE, 1u);
        expect(eq(one_range[0].front(), 0u));
        expect(eq(one_range[0].back(),  SIZE - 1));
    };

    "ranges_2"_test = [] {
        const auto two_ranges = init_uniform_ranges(SIZE, 2u);
        expect(eq(two_ranges[0].front(), 0u));
        expect(eq(two_ranges[0].back(),  0u));
        expect(eq(two_ranges[1].front(), 1u));
        expect(eq(two_ranges[1].back(),  SIZE - 1));
    };

    for(const size_t count : std::ranges::iota_view{3u, COUNTS}) {
        test("ranges_count_" + std::to_string(count)) = [count] {
            const auto more_ranges = init_uniform_ranges(SIZE, count);
            for(const size_t i : std::ranges::iota_view{0u, SIZE}) {
                expect(eq(more_ranges[i].front(), i));
                expect(eq(more_ranges[i].back(),  i));
            }
            for(const size_t i : std::ranges::iota_view{SIZE, more_ranges.size()}) {
                expect(eq(more_ranges[i].front(), SIZE));
                expect(eq(more_ranges[i].back(),  SIZE - 1));
            }
        };
    }
};

}
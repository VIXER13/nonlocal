#include "uniform_ranges.hpp"

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"uniform_range"> _ = [] {
    using namespace boost::ut;
    using namespace parallel;

    static constexpr size_t Size = 3;
    static constexpr size_t Counts = 6;
    
    "ranges_count_0"_test = [] {
        expect(throws([] { uniform_ranges(Size, 0u); })) <<
            "When creating 0 ranges, an exception should be thrown.";
    };

    for(const size_t count : std::ranges::iota_view{1u, Counts})
        test("ranges_count_" + std::to_string(count)) = [count] {
            const auto ranges = uniform_ranges(Size, count);
            expect(eq(ranges.size(), count)) << "Unexpected ranges number.";
        };

    "range_1"_test = [] {
        const auto one_range = uniform_ranges(Size, 1u);
        expect(eq(one_range[0].front(), 0u));
        expect(eq(one_range[0].back(),  Size - 1));
    };

    "ranges_2"_test = [] {
        const auto two_ranges = uniform_ranges(Size, 2u);
        expect(eq(two_ranges[0].front(), 0u));
        expect(eq(two_ranges[0].back(),  0u));
        expect(eq(two_ranges[1].front(), 1u));
        expect(eq(two_ranges[1].back(),  Size - 1));
    };

    for(const size_t count : std::ranges::iota_view{3u, Counts}) {
        test("ranges_" + std::to_string(count)) = [count] {
            const auto more_ranges = uniform_ranges(Size, count);
            for(const size_t i : std::ranges::iota_view{0u, Size}) {
                expect(eq(more_ranges[i].front(), i));
                expect(eq(more_ranges[i].back(),  i));
            }
            for(const size_t i : std::ranges::iota_view{Size, more_ranges.size()}) {
                expect(eq(more_ranges[i].front(), Size));
                expect(eq(more_ranges[i].back(),  Size - 1));
            }
        };
    }
};

}
#include <boost/ut.hpp>
#include <logger/logger.hpp>

using namespace logger;
using namespace boost::ut;

namespace
{
const suite<"logger"> _ = [] {
    "parse_level"_test = [] {
        // Test valid levels (case insensitive)
        expect(parse_level("off") == level::Off);
        expect(parse_level("OFF") == level::Off);
        expect(parse_level("error") == level::Error);
        expect(parse_level("ERROR") == level::Error);
        expect(parse_level("warning") == level::Warning);
        expect(parse_level("WARNING") == level::Warning);
        expect(parse_level("info") == level::Info);
        expect(parse_level("INFO") == level::Info);
        expect(parse_level("debug") == level::Debug);
        expect(parse_level("DEBUG") == level::Debug);
        expect(parse_level("trace") == level::Trace);
        expect(parse_level("TRACE") == level::Trace);

        // Test invalid levels
        expect(parse_level("invalid") == std::nullopt);
        expect(parse_level("") == std::nullopt);
    };

    "off"_test = [] {
        std::stringstream ss;
        auto& log = get(level::Off, std::make_unique<logger::stream_base>(ss));
        log << "This should not be logged.";
        expect(ss.str().empty());
        ss.clear();
    };
};
}
#include <metamath/functions/power.hpp>

#include <boost/ut.hpp>

#include <array>
#include <ranges>

namespace {

const boost::ut::suite<"power"> _ = [] {
    using namespace boost::ut;
    using namespace metamath::functions;

    static constexpr double Epsilon = std::numeric_limits<double>::epsilon();
    static constexpr auto Expected = std::array{0.03125, 0.0625, 0.125, 0.25, 0.5, 1., 2., 4., 8., 16., 32.};

    "static_power"_test = [] {
        expect(approx(power<-5>(2.), Expected[0],  Epsilon));
        expect(approx(power<-4>(2.), Expected[1],  Epsilon));
        expect(approx(power<-3>(2.), Expected[2],  Epsilon));
        expect(approx(power<-2>(2.), Expected[3],  Epsilon));
        expect(approx(power<-1>(2.), Expected[4],  Epsilon));
        expect(approx(power< 0>(2.), Expected[5],  Epsilon));
        expect(approx(power< 1>(2.), Expected[6],  Epsilon));
        expect(approx(power< 2>(2.), Expected[7],  Epsilon));
        expect(approx(power< 3>(2.), Expected[8],  Epsilon));
        expect(approx(power< 4>(2.), Expected[9],  Epsilon));
        expect(approx(power< 5>(2.), Expected[10], Epsilon));
    };

    "power_integer_exponent"_test = [] {
        for(const int exp : std::ranges::iota_view{-5, 6})
            expect(approx(power(2., exp), Expected[exp + 5], Epsilon));
    };

    "power_floating_point_exponent"_test = [] {
        for(const int exp : std::ranges::iota_view{-5, 6})
            expect(approx(power(2., double(exp)), Expected[exp + 5], Epsilon));
    };
};

}
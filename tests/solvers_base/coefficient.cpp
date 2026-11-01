#include <solvers/base/coefficient.hpp>

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"coefficient"> _ = [] {
    using namespace boost::ut;
    using namespace nonlocal;

    "evaluate_1d"_test = [] {
        const point<double, 1> x{ 2. };
        static constexpr coefficient_t<double, 1> Constant = 5.;
        static constexpr auto Expected = std::array{ 5., 6., 13. };
        const coefficient_t<double, 1> spatial = spatial_dependency<double, 1>{
            [](const point<double, 1>& x) constexpr noexcept { return 3. * x; }
        };
        const coefficient_t<double, 1> dependent = solution_dependency<double, 1>{
            [](const point<double, 1>& x, const double u) constexpr noexcept { return 3. * x + u; }
        };

        expect(eq(evaluate(Constant, x, 7.), 5.));
        expect(eq(evaluate(spatial, x, 7.), 6.));
        expect(eq(evaluate(dependent, x, 7.), 13.));
        expect(eq(evaluate(dependent, x, 9.), 15.));
        expect(evaluate(std::array{ Constant, spatial, dependent }, x, 7.) == Expected);
    };

    "evaluate_2d"_test = [] {
        const point<double, 2> x{ std::array{ 2., 3. } };
        static constexpr coefficient_t<double, 2> Constant = 5.;
        static constexpr auto Expected = std::array{ 5., 7., 14. };
        const coefficient_t<double, 2> spatial = spatial_dependency<double, 2>{
            [](const point<double, 2>& x) constexpr noexcept { return 2. * x[0] + x[1]; }
        };
        const coefficient_t<double, 2> dependent = solution_dependency<double, 2>{
            [](const point<double, 2>& x, const double u) constexpr noexcept { return 2. * x[0] + x[1] + u; }
        };

        expect(eq(evaluate(Constant, x, 7.), 5.));
        expect(eq(evaluate(spatial, x, 7.), 7.));
        expect(eq(evaluate(dependent, x, 7.), 14.));
        expect(eq(evaluate(dependent, x, 9.), 16.));
        expect(evaluate(std::array{ Constant, spatial, dependent }, x, 7.) == Expected);
    };
};

} // namespace

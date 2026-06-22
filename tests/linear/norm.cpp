#include <metamath/linear/norm.hpp>
#include <metamath/utils/constants.hpp>

#include <boost/ut.hpp>

namespace {

const boost::ut::suite<"norm"> _ = [] {
    using namespace boost::ut;
    using namespace metamath::linear;
    using T = double;

    static constexpr auto Epsilon = std::numeric_limits<T>::epsilon();
    static constexpr auto Inf = metamath::constants::Infinity<size_t>;

    "powered_norm_static_exponent"_test = []<class Container> {
        expect(approx(powered_norm<1>(Container{3.0, -4.0}), 7.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}), 25.0, Epsilon)); // By default Exp = 2
        expect(approx(powered_norm<3>(Container{3.0, -4.0}), 91.0, Epsilon));
        expect(approx(powered_norm<4>(Container{3.0, -4.0}), 337.0, Epsilon));
        expect(eq(powered_norm<Inf>(Container{3.0, -4.0}), std::numeric_limits<T>::infinity()));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};

    "powered_norm_integer_exponent"_test = []<class Container> {
        expect(approx(powered_norm(Container{3.0, -4.0}, 1), 7.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 2), 25.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 3), 91.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 4), 337.0, Epsilon));
        expect(eq(powered_norm(Container{3.0, -4.0}, Inf), std::numeric_limits<T>::infinity()));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};

    "powered_norm_floating_point_exponent"_test = []<class Container> {
        expect(approx(powered_norm(Container{3.0, -4.0}, 1.0), 7.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 2.0), 25.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 3.0), 91.0, Epsilon));
        expect(approx(powered_norm(Container{3.0, -4.0}, 4.0), 337.0, Epsilon));
        expect(eq(powered_norm(Container{3.0, -4.0}, T(Inf)), std::numeric_limits<T>::infinity()));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};

    "norm_static_exponent"_test = []<class Container> {
        expect(approx(norm<1>(Container{3.0, -4.0}), 7.0, Epsilon));
        expect(approx(norm(Container{3.0, -4.0}), 5.0, Epsilon)); // By default Exp = 2
        expect(approx(norm<3>(Container{3.0, -4.0}), std::cbrt(91.0), Epsilon));
        expect(approx(norm<4>(Container{3.0, -4.0}), std::sqrt(std::sqrt(337.0)), Epsilon));
        expect(eq(norm<Inf>(Container{3.0, -4.0}), 4.0));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};

    "norm_integer_exponent"_test = []<class Container> {
        expect(approx(norm(Container{3.0, -4.0}), 5.0, Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 1), 7.0, Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 3), std::cbrt(91.0), Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 4), std::sqrt(std::sqrt(337.0)), Epsilon));
        expect(eq(norm(Container{3.0, -4.0}, Inf), std::numeric_limits<T>::infinity()));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};

    "norm_floating_point_exponent"_test = []<class Container> {
        expect(approx(norm(Container{3.0, -4.0}), 5.0, Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 1.0), 7.0, Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 3.0), std::cbrt(91.0), Epsilon));
        expect(approx(norm(Container{3.0, -4.0}, 4.0), std::sqrt(std::sqrt(337.0)), Epsilon));
        expect(eq(norm(Container{3.0, -4.0}, T(Inf)), std::numeric_limits<T>::infinity()));
    } | std::tuple<std::array<T, 2>, std::vector<T>>{};
};

}
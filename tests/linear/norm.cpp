#include <metamath/linear/norm.hpp>
#include <metamath/utils/constants.hpp>

#include <boost/ut.hpp>

namespace {

using namespace boost::ut;
using namespace metamath::linear;
using T = double;
using containers = std::tuple<std::array<T, 2>, std::vector<T>, std::array<std::array<T, 2>, 2>, std::vector<std::array<T, 2>>>;

constexpr auto Epsilon = std::numeric_limits<T>::epsilon();
constexpr auto Inf = metamath::constants::Infinity<size_t>;

template<class Container>
Container init_container() {
    if constexpr (metamath::types::is_array_v<std::ranges::range_value_t<Container>>)
        return {std::array<T, 2>{3.0, -4.0}, std::array<T, 2>{5.0, -6.0}};
    else
        return {3.0, -4.0};
}

template<class Container>
constexpr std::array<T, 5> excpected_powered_norms() noexcept {
    return metamath::types::is_array_v<std::ranges::range_value_t<Container>>
         ? std::array{18.0, 86.0, 432.0, 2258.0, std::numeric_limits<T>::infinity()}
         : std::array{7.0, 25.0, 91.0, 337.0, std::numeric_limits<T>::infinity()};
}

template<class Container>
std::array<T, 5> excpected_norms() noexcept {
    return metamath::types::is_array_v<std::ranges::range_value_t<Container>>
         ? std::array{18.0, std::sqrt(86.0), std::cbrt(432.0), std::sqrt(std::sqrt(2258.0)), std::numeric_limits<T>::infinity()}
         : std::array{7.0, 5.0, std::cbrt(91.0), std::sqrt(std::sqrt(337.0)), std::numeric_limits<T>::infinity()};
}

const suite<"norm"> _ = [] {
    "powered_norm_static_exponent"_test = []<class Container> {
        constexpr auto Expected = excpected_powered_norms<Container>();
        expect(approx(powered_norm<1>(init_container<Container>()), Expected[0], Epsilon));
        expect(approx(powered_norm(init_container<Container>()), Expected[1], Epsilon)); // By default Exp = 2
        expect(approx(powered_norm<3>(init_container<Container>()), Expected[2], Epsilon));
        expect(approx(powered_norm<4>(init_container<Container>()), Expected[3], Epsilon));
        expect(eq(powered_norm<Inf>(init_container<Container>()), Expected[4]));
    } | containers{};

    "powered_norm_integer_exponent"_test = []<class Container> {
        constexpr auto Expected = excpected_powered_norms<Container>();
        expect(approx(powered_norm(init_container<Container>(), 1), Expected[0], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 2), Expected[1], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 3), Expected[2], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 4), Expected[3], Epsilon));
        expect(eq(powered_norm(init_container<Container>(), Inf), Expected[4]));
    } | containers{};

    "powered_norm_floating_point_exponent"_test = []<class Container> {
        constexpr auto Expected = excpected_powered_norms<Container>();
        expect(approx(powered_norm(init_container<Container>(), 1.0), Expected[0], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 2.0), Expected[1], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 3.0), Expected[2], Epsilon));
        expect(approx(powered_norm(init_container<Container>(), 4.0), Expected[3], Epsilon));
        expect(eq(powered_norm(init_container<Container>(), T(Inf)), Expected[4]));
    } | containers{};

    "norm_static_exponent"_test = []<class Container> {
        using RT = std::ranges::range_value_t<Container>;
        const auto expected = excpected_norms<Container>();
        expect(approx(norm<1>(init_container<Container>()), expected[0], Epsilon));
        expect(approx(norm(init_container<Container>()), expected[1], Epsilon)); // By default Exp = 2
        expect(approx(norm<3>(init_container<Container>()), expected[2], 9e-16));
        expect(approx(norm<4>(init_container<Container>()), expected[3], Epsilon));
        expect(eq(norm<Inf>(init_container<Container>()), std::is_same_v<RT, T> ? 4.0 : 6.0));
    } | containers{};

    "norm_integer_exponent"_test = []<class Container> {
        const auto expected = excpected_norms<Container>();
        expect(approx(norm(init_container<Container>(), 1), expected[0], Epsilon));
        expect(approx(norm(init_container<Container>(), 2), expected[1], Epsilon));
        expect(approx(norm(init_container<Container>(), 3), expected[2], 2e-15));
        expect(approx(norm(init_container<Container>(), 4), expected[3], Epsilon));
        expect(eq(norm(init_container<Container>(), Inf), expected[4]));
    } | containers{};

    "norm_floating_point_exponent"_test = []<class Container> {
        const auto expected = excpected_norms<Container>();
        expect(approx(norm(init_container<Container>(), 1.0), expected[0], Epsilon));
        expect(approx(norm(init_container<Container>(), 2.0), expected[1], Epsilon));
        expect(approx(norm(init_container<Container>(), 3.0), expected[2], 2e-15));
        expect(approx(norm(init_container<Container>(), 4.0), expected[3], Epsilon));
        expect(eq(norm(init_container<Container>(), T(Inf)), expected[4]));
    } | containers{};
};

}
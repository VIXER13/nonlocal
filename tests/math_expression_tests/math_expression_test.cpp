#include <math_expression/math_expression.hpp>

#include <boost/ut.hpp>

#include <numbers>
#include <limits>

namespace {

using namespace boost::ut;
using namespace formula;
using namespace std::numbers;
using namespace std::string_literals;
using T = double;

const suite<"formula"> _ = [] {
    "polish_notation_general"_test = [] {
        auto test = math_expression<T>{"x : -x"};
        expect(eq(test.to_polish(), "x ~"s));

        test = math_expression<T>{"x:-x"};
        expect(eq(test.to_polish(), "x ~"s));

        test = math_expression<T>{"x : x + x - x / x * x"};
        expect(eq(test.to_polish(), "x x + x x / x * -"s));

        test = math_expression<T>{"x : sin(x)"};
        expect(eq(test.to_polish(), "x sin"s));
        expect(lt(std::abs(test({ pi / 4 }) - std::sin(pi / 4)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : cos(x)"};
        expect(eq(test.to_polish(), "x cos"s));
        expect(lt(std::abs(test({ pi / 4 }) - std::cos(pi / 4)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : tan(x)"};
        expect(eq(test.to_polish(), "x tan"s));
        expect(lt(std::abs(test({ pi / 4 }) - std::tan(pi / 4)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : atan(x)"};
        expect(eq(test.to_polish(), "x atan"s));
        expect(lt(std::abs(test({ pi / 4 }) - std::atan(pi / 4)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : abs(x)"};
        expect(eq(test.to_polish(), "x abs"s));
        expect(lt(std::abs(test({ -pi / 4 }) - std::abs(-pi / 4)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : sign(x)"};
        expect(eq(test.to_polish(), "x sign"s));
        expect(lt(std::abs(test({ -pi / 4 }) + 1), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : sqr(x)"};
        expect(eq(test.to_polish(), "x sqr"s));
        expect(lt(std::abs(test({ std::sqrt(2) }) - 2), std::numeric_limits<double>::epsilon() * 3));

        test = math_expression<T>{"x : sqrt(x)"};
        expect(eq(test.to_polish(), "x sqrt"s));
        expect(lt(std::abs(test({ 5. }) - std::sqrt(5)), std::numeric_limits<double>::epsilon() * 3));

        test = math_expression<T>{"x : log(x)"};
        expect(eq(test.to_polish(), "x log"s));
        expect(lt(std::abs(test({ e }) - 1), std::numeric_limits<double>::epsilon()));

        // remark: Euler Gamma function
        // for each natural number i = 1, 2, 3, ... tgamma(i) = (i - 1)! 
        test = math_expression<T>{"x : tgamma(x)"};
        expect(eq(test.to_polish(), "x tgamma"s));
        expect(lt(std::abs(test({ 6. }) - 120.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : lgamma(x)"};
        expect(eq(test.to_polish(), "x lgamma"s));
        expect(lt(std::abs(test({ 124.23 }) - std::lgamma(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : exp2(x)"};
        expect(eq(test.to_polish(), "x exp2"s));
        expect(lt(std::abs(test({ 4. }) - 16.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : expm1(x)"};
        expect(eq(test.to_polish(), "x expm1"s));
        expect(lt(std::abs(test({ 0. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : log10(x)"};
        expect(eq(test.to_polish(), "x log10"s));
        expect(lt(std::abs(test({ 100. }) - 2.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : log2(x)"};
        expect(eq(test.to_polish(), "x log2"s));
        expect(lt(std::abs(test({ 4. }) - 2.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : log1p(x)"};
        expect(eq(test.to_polish(), "x log1p"s));
        expect(lt(std::abs(test({ e - 1. }) - 1.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : cbrt(x)"};
        expect(eq(test.to_polish(), "x cbrt"s));
        expect(approx(test({ T{27} }), std::cbrt(T{27}), 3 * std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : asin(x)"};
        expect(eq(test.to_polish(), "x asin"s));
        expect(lt(std::abs(test({ 1. }) - pi / 2.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : acos(x)"};
        expect(eq(test.to_polish(), "x acos"s));
        expect(lt(std::abs(test({ 0. }) - pi / 2.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : sinh(x)"};
        expect(eq(test.to_polish(), "x sinh"s));
        expect(lt(std::abs(test({ 0. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : cosh(x)"};
        expect(eq(test.to_polish(), "x cosh"s));
        expect(lt(std::abs(test({ 0. }) - 1.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : tanh(x)"};
        expect(eq(test.to_polish(), "x tanh"s));
        expect(lt(std::abs(test({ 0. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : asinh(x)"};
        expect(eq(test.to_polish(), "x asinh"s));
        expect(lt(std::abs(test({ 0. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : acosh(x)"};
        expect(eq(test.to_polish(), "x acosh"s));
        expect(lt(std::abs(test({ 1. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : atanh(x)"};
        expect(eq(test.to_polish(), "x atanh"s));
        expect(lt(std::abs(test({ 0. }) - 0.), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : erf(x)"};
        expect(eq(test.to_polish(), "x erf"s));
        expect(lt(std::abs(test({ 124.23 }) - std::erf(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : erfc(x)"};
        expect(eq(test.to_polish(), "x erfc"s));
        expect(lt(std::abs(test({ 124.23 }) - std::erfc(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : trunc(x)"};
        expect(eq(test.to_polish(), "x trunc"s));
        expect(lt(std::abs(test({ 124.23 }) - std::trunc(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : floor(x)"};
        expect(eq(test.to_polish(), "x floor"s));
        expect(lt(std::abs(test({ 124.23 }) - std::floor(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : ceil(x)"};
        expect(eq(test.to_polish(), "x ceil"s));
        expect(lt(std::abs(test({ 124.23 }) - std::ceil(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : round(x)"};
        expect(eq(test.to_polish(), "x round"s));
        expect(lt(std::abs(test({ 124.23 }) - std::round(124.23)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : x^(0.3415234)"};
        expect(eq(test.to_polish(), "x 0.341523 ^"s)); // only 6 signs after comma!
        expect(lt(std::abs(test({ e }) - std::pow(e, 0.3415234)), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x : 0.411313 + .5 - x"};
        expect(eq(test.to_polish(), "0.411313 0.500000 + x -"s));
        expect(lt(std::abs(test({ 0.411313 }) - 0.5), std::numeric_limits<double>::epsilon()));

        // trivial literal arithmetics
        test = math_expression<T>{"x : -5 + 56.23424 - .51241 / 0.4321 * 4 / 10"};
        expect(eq(test.to_polish(), "5.000000 ~ 56.234240 + 0.512410 0.432100 / 4.000000 * 10.000000 / -"s));
        expect(lt(std::abs(test({ 100. }) - 50.7599), 1e-5));

        // general routine, examples
        test = math_expression<T>{"v1 v2 v3 p1 p2 p3 m : (p1^2 + p2^2 + p3^2) / (2 * m) + .5 * (v1^2 + v2^2 + v3^2) "};
        expect(eq(test.variables_count(), 7));
        expect(eq(test.to_polish(), "p1 2.000000 ^ p2 2.000000 ^ + p3 2.000000 ^ + 2.000000 m * / 0.500000 v1 2.000000 ^ v2 2.000000 ^ + v3 2.000000 ^ + * +"s));
        expect(lt(std::abs(test({ 1., 1., 1., 2., 2., 2., 2. }) - 4.5), 1e-8));

        test = math_expression<T>{"x a b : (x - 1)^(a - 1) * (x + 1)^(b + 1)"};
        expect(eq(test.variables_count(), 3));
        expect(eq(test.to_polish(), "x 1.000000 - a 1.000000 - ^ x 1.000000 + b 1.000000 + ^ *"s));
        expect(lt(std::abs(test({ 1., 1., 1. }) - 4), std::numeric_limits<double>::epsilon()));
        expect(lt(std::abs(test({ 1., 100., 100. })), std::numeric_limits<double>::epsilon()));
        expect(lt(std::abs(test({ -1., 100., 100. })), std::numeric_limits<double>::epsilon()));
        expect(lt(std::abs(test({ 3., 3., 2. }) - 256), std::numeric_limits<double>::epsilon()));

        test = math_expression<T>{"x y t : x * cos(y) / exp(t) + 10"};
        expect(eq(test.variables_count(), 3));
        expect(eq(test.to_polish(), "x y cos * t exp / 10.000000 +"s));
        expect(std::abs(test({ 2., pi, 10. }) - 9.99991) < 1e-6);

        test = math_expression<T>{"x y : 400.0 * exp((x - 1.0)^2 * -400.0)"};
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "400.000000 x 1.000000 - 2.000000 ^ 400.000000 ~ * exp *"s));
        expect(std::abs(test({ 1., 100. }) - 400.) < 1e-6);
        expect(std::abs(test({ 1.13, 100. }) - 0.463692) < 1e-6);
    };

    "polish_notation_throws"_test = [] {
        const auto& operator_priority = utils::get_operator_priority();

        // Wrong variables format. Symbol ':' is required after variables initialization.
        expect(throws([]() { math_expression<T>{"x y z x * y * z"}; }));

        // Wrong variables format. Variables must start with latin letter, variable cannot start with a number
        expect(throws([]() { math_expression<T>{"1x 2x 3x: 1x + 2x + 3x)"}; }));
        expect(nothrow([]() { math_expression<T>{"x1 x2 x3: x1 + x2 + x3"}; }));

        // Wrong expression format. Formula after ':' is required.
        expect(throws([]() { math_expression<T>{"x y z: "}; }));

        // Everything alright. Variables are not required if not used.
        expect(nothrow([]() { math_expression<T>{": 2 * 2"}; }));
        expect(nothrow([]() { math_expression<T>{"x : 2 * 2"}; }));

        // TODO : implement verification for this case, when variable in infix notation is not used 
        auto test = math_expression<T>{"x : 2 * 2 * y"};
        //expect(test.to_polish() == "22**");
        // Must be throws
        expect(nothrow([]() { math_expression<T>{"x : 2 * 2 * y"}; }));
        // Must be throws, no functionality related with symbol "Г"
        expect(nothrow([]() { math_expression<T>{"x y z: Г(x) * Г(y) * Г(z)"}; }));

        // Invalid variable designation. Variable name <" + variable + "> is unavailable.
        for (const auto& node : operator_priority) 
            expect(throws([op = node.first]() { math_expression<T>{op + " : 2 * 2"}; }));
        
        // Wrong expression format. The expression contains open parentheses
        expect(throws([]() { math_expression<T>{"x y z: exp(x*y*z"}; }));
        expect(throws([]() { math_expression<T>{"x y z: exp x*y*z)"}; }));

        // Wrong expression format. The expression can not be started or finished with dot. Dots are only allowed in number representation.
        // TODO : 0. format
        expect(throws([]() { math_expression<T>{"x y z: .5"}; }));
        expect(throws([]() { math_expression<T>{"x y z: 0."}; }));
        // must be nothrow
        expect(throws([]() { math_expression<T>{"x y z: 0. * x / y"}; }));
        expect(nothrow([]() { math_expression<T>{"x y z: 0.5"}; }));
        expect(nothrow([]() { math_expression<T>{"x y z: 0.0"}; }));
        expect(nothrow([]() { math_expression<T>{"x y z: x* .5 / y"}; }));
    };
};

}
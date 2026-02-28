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
constexpr auto eps = std::numeric_limits<T>::epsilon();
const suite<"formula"> _ = [] {
    "math_expression_literal"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization

        expect(nothrow([&test] {test = math_expression<T>{" : 1"};})) << fatal;
        expect(eq(test.variables_count(), 0));
        expect(eq(test.to_polish(), "1.000000"s));
        expect(lt(std::abs(test({}) - 1.0), eps));

        expect(nothrow([&test] {test = math_expression<T>{" : .1"};})) << fatal;
        expect(eq(test.variables_count(), 0));
        expect(eq(test.to_polish(), "0.100000"s));
        expect(lt(std::abs(test({}) - 0.1), eps));

        expect(nothrow([&test] {test = math_expression<T>{" : 1."};})) << fatal;
        expect(eq(test.variables_count(), 0));
        expect(eq(test.to_polish(), "1.000000"s));
        expect(lt(std::abs(test({}) - 1.0), eps));

        expect(nothrow([&test] {test = math_expression<T>{" : 0.123456789"};})) << fatal;
        expect(eq(test.variables_count(), 0));
        expect(eq(test.to_polish(), "0.123457"s));
        expect(lt(std::abs(test({}) - 0.123456789), eps));
    };
    "math_expression_unary"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization
        for( auto& [name, func] : math_expression<T>::unary_operators()) {
            expect(nothrow([&test, &name] {test = math_expression<T>{"x : " + name + "(x)"};})) << fatal;
            expect(eq(test.to_polish(), "x " + name)) << name << " operator is not correctly parsed.";
            for(auto x : {0.1, 0.2, pi/4, 1., 2., 5.})
                if(!std::isnan(func(x)) && !std::isinf(func(x)))
                    expect(lt(std::abs(test({x}) - func(x)), eps)) << name << "operator is not correctly computed at x =" << x;
        }
    };
    "math_expression_binary"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization
        for( auto& [name, func] : math_expression<T>::binary_operators()) {
            expect(nothrow([&test, &name] {test = math_expression<T>{"x y : x " + name + " y"};})) << fatal;
            expect(eq(test.to_polish(), "x y " + name)) << name << " operator is not correctly parsed.";
            for(auto x : {0.1, 0.2, pi/4, 1., 2., 5.})
                for(auto y : {0.1, 0.2, pi/4, 1., 2., 5.})
                if(!std::isnan(func(x, y)) && !std::isinf(func(x, y)))
                    expect(lt(std::abs(test({x, y}) - func(x, y)), eps)) << name << "operator is not correctly computed at x =" << x << ", y =" << y;
        }
    };
    "math_expression_complex"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization

        expect(nothrow([&test] {test = math_expression<T>{"x : 0.411313 + .5 - x"};})) << fatal;
        expect(eq(test.variables_count(), 1));
        expect(eq(test.to_polish(), "0.411313 0.500000 + x -"s));
        expect(lt(std::abs(test({ 0.411313 }) - 0.5), eps));
        
        expect(nothrow([&test] {test = math_expression<T>{"x : -5 + 56.23424 - .51241 / 0.4321 * 4 / 10"};})) << fatal;
        expect(eq(test.variables_count(), 1));
        expect(eq(test.to_polish(), "5.000000 ~ 56.234240 + 0.512410 0.432100 / 4.000000 * 10.000000 / -"s));
        expect(lt(std::abs(test({ 100. }) - 50.7599), 1e-5));

        expect(nothrow([&test] {test = math_expression<T>{"v1 v2 v3 p1 p2 p3 m : (p1^2 + p2^2 + p3^2) / (2 * m) + .5 * (v1^2 + v2^2 + v3^2) "};})) << fatal;
        expect(eq(test.variables_count(), 7));
        expect(eq(test.to_polish(), "p1 2.000000 ^ p2 2.000000 ^ + p3 2.000000 ^ + 2.000000 m * / 0.500000 v1 2.000000 ^ v2 2.000000 ^ + v3 2.000000 ^ + * +"s));
        expect(lt(std::abs(test({ 1., 1., 1., 2., 2., 2., 2. }) - 4.5), eps));

        expect(nothrow([&test] {test = math_expression<T>{"x a b : (x - 1)^(a - 1) * (x + 1)^(b + 1)"};})) << fatal;
        expect(eq(test.variables_count(), 3));
        expect(eq(test.to_polish(), "x 1.000000 - a 1.000000 - ^ x 1.000000 + b 1.000000 + ^ *"s));
        expect(lt(std::abs(test({ 1., 1., 1. }) - 4), eps));
        expect(lt(std::abs(test({ 1., 100., 100. })), eps));
        expect(lt(std::abs(test({ -1., 100., 100. })), eps));
        expect(lt(std::abs(test({ 3., 3., 2. }) - 256), eps));

        expect(nothrow([&test] {test = math_expression<T>{"x y t : x * cos(y) / exp(t) + 10"};})) << fatal;
        expect(eq(test.variables_count(), 3));
        expect(eq(test.to_polish(), "x y cos * t exp / 10.000000 +"s));
        expect(std::abs(test({ 2., pi, 10. }) - 9.99991) < 1e-6);

        expect(nothrow([&test] {test = math_expression<T>{"x y : 400.0 * exp((x - 1.0)^2 * -400.0)"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "400.000000 x 1.000000 - 2.000000 ^ 400.000000 ~ * exp *"s));
        expect(std::abs(test({ 1., 100. }) - 400.) < 1e-6);
        expect(std::abs(test({ 1.13, 100. }) - 0.463692) < 1e-6);

        expect(nothrow([&test] {test = math_expression<T>{"x t : exp(-t) * (sin(x)^2 + cos(x)^2)"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "t ~ exp x sin 2.000000 ^ x cos 2.000000 ^ + *"s));
        expect(lt(std::abs(test({ 0.3, 1.7 }) - std::exp(-1.7)), eps));

        expect(nothrow([&test] {test = math_expression<T>{"x y : sqrt(x^2 + y^2) + abs(x - y) + (-sign(x*y))"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "x 2.000000 ^ y 2.000000 ^ + sqrt x y - abs + x y * sign ~ +"s));
        expect(lt(std::abs(test({ 3., 4. }) - 5.), eps));

        expect(nothrow([&test] {test = math_expression<T>{"x y : atan2(y, x) + log1p(abs(x)) * expm1(y)"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "y x atan2 x abs log1p y expm1 * +"s));
        expect(lt(std::abs(test({ 1., 0.5 }) - (std::atan2(0.5, 1.0) + std::log1p(std::abs(1.0)) * std::expm1(0.5))), eps));

        // regression: binary '-' after function call (used to be mis-parsed as unary)
        expect(nothrow([&test] {test = math_expression<T>{"x y : abs(x) - y"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "x abs y -"s));
        expect(lt(std::abs(test({ -3., 4. }) - (std::abs(-3.) - 4.)), eps));

        // regression: comma must flush operators inside function arguments
        expect(nothrow([&test] {test = math_expression<T>{"x y : atan2(y + 1, x * 2)"};})) << fatal;
        expect(eq(test.variables_count(), 2));
        expect(eq(test.to_polish(), "y 1.000000 + x 2.000000 * atan2"s));
        expect(lt(std::abs(test({ 3., 4. }) - std::atan2(4. + 1., 3. * 2.)), eps));

        expect(nothrow([&test] {test = math_expression<T>{"u v w : (u + v - w)^3 / (1 + sqr(u - v))"};})) << fatal;
        expect(eq(test.variables_count(), 3));
        expect(eq(test.to_polish(), "u v + w - 3.000000 ^ 1.000000 u v - sqr + /"s));
        expect(lt(std::abs(test({ 2., 0.5, 1. }) - (std::pow(2.0 + 0.5 - 1.0, 3.0) / (1.0 + (2.0 - 0.5) * (2.0 - 0.5)))), eps));

        expect(nothrow([&test] {test = math_expression<T>{"x : exp2(log2(x)) + expm1(log1p(x))"};})) << fatal;
        expect(eq(test.variables_count(), 1));
        expect(eq(test.to_polish(), "x log2 exp2 x log1p expm1 +"s));
        expect(lt(std::abs(test({ 0.25 }) - 0.5), eps));
    };

    "math_expression_throws"_test = [] {
        const auto& operator_priority = utils::get_operator_priority();

        // Wrong variables format. Symbol ':' is required after variables initialization.
        expect(throws([]() { math_expression<T>{"x y z x * y * z"}; }));

        // Wrong variables format. Variables must start with latin letter, variable cannot start with a number
        expect(throws([]() { math_expression<T>{"1x 2x 3x: 1x + 2x + 3x)"}; }));

        // Wrong expression format. Formula after ':' is required.
        expect(throws([]() { math_expression<T>{"x y z: "}; }));

        // Unknown variable used in expression.
        expect(throws([]() { math_expression<T>{"x : 2 * 2 * y"}; }));

        // Wrong expression format. non-ANSI "Г" present.
        expect(throws([]() { math_expression<T>{"x y z: Г(x) * Г(y) * Г(z)"}; }));

        // Invalid variable designation. Variable name <" + variable + "> is unavailable.
        for (const auto& node : operator_priority) 
            expect(throws([op = node.first]() { math_expression<T>{op + " : 2 * 2"}; }));
        
        // Wrong expression format. The expression contains open parentheses
        expect(throws([]() { math_expression<T>{"x y z: exp(x*y*z"}; }));
        expect(throws([]() { math_expression<T>{"x y z: exp x*y*z)"}; }));

        // Wrong variables format. Variables must be declared as symbols before ':' separator.
        expect(throws([]() { math_expression<T>{"x 1 : x"}; }));
        expect(throws([]() { math_expression<T>{"( x ) : x"}; }));

        // Tokenizer edge cases
        // Invalid numeric construction: digit sequence immediately followed by symbol characters.
        expect(throws([]() { math_expression<T>{"x : 1x + 2"}; }));
        // Invalid dots in number representation.
        expect(throws([]() { math_expression<T>{"x : 1..2"}; }));
        // Unexpected / unsupported character can lead to an empty-token tokenizer failure.
        expect(throws([]() { math_expression<T>{"x : x $ 1"}; }));

        // Wrong expression format. Unexpected separator in expression part.
        expect(throws([]() { math_expression<T>{"x : x : 1"}; }));
        expect(throws([]() { math_expression<T>{"x y : atan2(x|y, 1)"}; }));

        // Unknown function/operator in expression.
        expect(throws([]() { math_expression<T>{"x : foo(x)"}; }));

        // Wrong number of variables at evaluation time.
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization
        expect(nothrow([&test] {test =  math_expression<T>{"x y : x + y"};})) << fatal;
        expect(throws([&test]() { (void)test({1.}); }));
        expect(throws([&test]() { (void)test({1., 2., 3.}); }));
    };
};

}
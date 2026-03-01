#include <math_expression/math_expression.hpp>
#include <math_expression/utils.hpp>

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

// Some ugly code to simplify testing
using std::abs; using std::acos; using std::acosh; using std::asin; using std::asinh; 
using std::atan; using std::atanh; using std::cbrt; using std::ceil; using std::cos; 
using std::cosh; using std::erf; using std::erfc; using std::exp; using std::exp2; 
using std::expm1; using std::floor; using std::lgamma; using std::log; using std::log10; 
using std::log1p; using std::log2; using std::round; using std::sin; using std::sinh; 
using std::sqrt; using std::tan; using std::tanh; using std::tgamma; using std::trunc; 
using std::pow; using std::hypot; using std::fmod; using std::min; using std::max;
static constexpr auto sqr = [](const T value) { return value * value; };
static constexpr auto sign = [](const T value) noexcept { return T((value > 0) - (value < 0)); };
#define _COUNT(...) (std::ranges::count(#__VA_ARGS__, ',') + 1)
#define _FUNC(exp,...) \
[](std::array<T, _COUNT(__VA_ARGS__)> arr) -> T { \
    const auto& [__VA_ARGS__] = arr;\
    if constexpr (std::ranges::count(#exp, '^') > 0) return std::nan("1"); \
    else return exp;\
}
#define PRESETUP_EXPR(exp,...) \
        auto repr = #exp; \
        expect(nothrow([&test, &repr] {test = math_expression<T>{#__VA_ARGS__ ":" #exp};})) << repr << "expression ctor throws" << fatal; \
        expect(eq(test.variables_count(), _COUNT(__VA_ARGS__))) << repr << "expression has wrong variables count"
#define SETUP_EXPR(exp,...) \
        PRESETUP_EXPR(exp, __VA_ARGS__); \
        auto func = _FUNC(exp, __VA_ARGS__)
#define CHECK_VALUE(...) \
        expect(lt(std::abs(test({__VA_ARGS__}) - func({__VA_ARGS__})), eps)) << repr << "expression is not correctly computed at values = {" #__VA_ARGS__ "}";
#define CHECK_NOTATION(expected) \
        expect(eq(test.to_polish(), expected ## s)) << repr << "expression is not correctly parsed"

const suite<"formula"> _ = [] {
    "tokenize"_test = [] {
        using formula::utils::token_t;
        using enum formula::utils::token_t::type_t;
        using formula::utils::tokenize;

        expect(eq(tokenize("").size(), 0));
        expect(eq(tokenize(" \t\n\r\v\f ").size(), 0));

        {
            const auto tokens = tokenize("x y : sin(x) + 2.5*y");
            expect(eq(tokens.size(), 11));
            expect(eq(tokens[0], token_t{Symbol, "x"}));
            expect(eq(tokens[1], token_t{Symbol, "y"}));
            expect(eq(tokens[2], token_t{Separator, ":"}));
            expect(eq(tokens[3], token_t{Symbol, "sin"}));
            expect(eq(tokens[4], token_t{ParenthesisLeft, "("}));
            expect(eq(tokens[5], token_t{Symbol, "x"}));
            expect(eq(tokens[6], token_t{ParenthesisRight, ")"}));
            expect(eq(tokens[7], token_t{Operator, "+"}));
            expect(eq(tokens[8], token_t{Number, "2.5"}));
            expect(eq(tokens[9], token_t{Operator, "*"}));
            expect(eq(tokens[10], token_t{Symbol, "y"}));
        }

        {
            const auto tokens = tokenize("atan2(x, y)");
            expect(eq(tokens.size(), 6));
            expect(eq(tokens[0], token_t{Symbol, "atan2"}));
            expect(eq(tokens[1], token_t{ParenthesisLeft, "("}));
            expect(eq(tokens[2], token_t{Symbol, "x"}));
            expect(eq(tokens[3], token_t{Separator, ","}));
            expect(eq(tokens[4], token_t{Symbol, "y"}));
            expect(eq(tokens[5], token_t{ParenthesisRight, ")"}));
        }

        {
            const auto tokens = tokenize(".1 1.");
            expect(eq(tokens.size(), 2));
            expect(eq(tokens[0], token_t{Number, ".1"}));
            expect(eq(tokens[1], token_t{Number, "1."}));
        }
    };

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
            expect(nothrow([&test, &name] {test = math_expression<T>{"x : " + name + "(x)"};})) << name << "operator ctor throws" << fatal;
            expect(eq(test.to_polish(), "x " + name)) << name << "operator is not correctly parsed.";
            for(auto x : {0.1, 0.2, pi/4, 1., 2., 5.})
                if(!std::isnan(func(x)) && !std::isinf(func(x)))
                    expect(lt(std::abs(test({x}) - func(x)), eps)) << name << "operator is not correctly computed at x =" << x;
        }
    };
    "math_expression_binary"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization
        for(auto& [name, func] : math_expression<T>::binary_operators()) {
            expect(nothrow([&test, &name] {test = math_expression<T>{"x y : x " + name + " y"};})) << name << "operator ctor throws" << fatal;
            expect(eq(test.to_polish(), "x y " + name)) << name << "operator is not correctly parsed.";
            for(auto x : {0.1, 0.2, pi/4, 1., 2., 5.})
                for(auto y : {0.1, 0.2, pi/4, 1., 2., 5.})
                if(!std::isnan(func(x, y)) && !std::isinf(func(x, y)))
                    expect(lt(std::abs(test({x, y}) - func(x, y)), eps)) << name << "operator is not correctly computed at x =" << x << ", y =" << y;

            if(name.size() == 1) // skip single symbol operators
                continue;
            expect(nothrow([&test, &name] {test = math_expression<T>{"x y : " + name + "(x, y)"};})) << name << "operator ctor throws" << fatal;
            expect(eq(test.to_polish(), "x y " + name)) << name << "operator is not correctly parsed.";
            for(auto x : {0.1, 0.2, pi/4, 1., 2., 5.})
                for(auto y : {0.1, 0.2, pi/4, 1., 2., 5.})
                if(!std::isnan(func(x, y)) && !std::isinf(func(x, y)))
                    expect(lt(std::abs(test({x, y}) - func(x, y)), eps)) << name << "operator is not correctly computed at x =" << x << ", y =" << y;

        }
    };
    "math_expression_complex"_test = [] {
        auto test = math_expression<T>{" : 0"}; // dummy expression for initialization
        {
            SETUP_EXPR(0.411313 + .5 - x, x);
            CHECK_NOTATION("0.411313 0.500000 + x -");
            CHECK_VALUE(0.411313);
        }

        {
            SETUP_EXPR(-5 + 56.23424 - .51241 / 0.4321 * 4 / 10, x);
            CHECK_NOTATION("5.000000 ~ 56.234240 + 0.512410 0.432100 / 4.000000 * 10.000000 / -");
            CHECK_VALUE(100);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR((p1^2 + p2^2 + p3^2) / (2 * m) + .5 * (v1^2 + v2^2 + v3^2), v1, v2, v3, p1, p2, p3, m);
            auto func = [](std::array<T, 7> arr) -> T {
                const auto& [v1, v2, v3, p1, p2, p3, m] = arr;
                return (sqr(p1) + sqr(p2) + sqr(p3)) / (2 * m) + .5 * (sqr(v1) + sqr(v2) + sqr(v3));
            };
            CHECK_NOTATION("p1 2.000000 ^ p2 2.000000 ^ + p3 2.000000 ^ + 2.000000 m * / 0.500000 v1 2.000000 ^ v2 2.000000 ^ + v3 2.000000 ^ + * +");
            CHECK_VALUE({ 1., 1., 1., 2., 2., 2., 2. });
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR((x - 1)^(a - 1) * (x + 1)^(b + 1), x, a, b);
            auto func = [](std::array<T, 3> arr) -> T {
                const auto& [x, a, b] = arr;
                return pow(x - 1, a - 1) * pow(x + 1, b + 1);
            };
            CHECK_NOTATION("x 1.000000 - a 1.000000 - ^ x 1.000000 + b 1.000000 + ^ *");
            CHECK_VALUE({ 1., 1., 1. });
            CHECK_VALUE({ 1., 100., 100. });
            CHECK_VALUE({ -1., 100., 100. });
            CHECK_VALUE({ 3., 3., 2. });
        }

        {
            SETUP_EXPR(x * cos(y) / exp(t) + 10, x, y, t);
            CHECK_NOTATION("x y cos * t exp / 10.000000 +");
            CHECK_VALUE(2., pi, 10.);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR(400.0 * exp((x - 1.0)^2 * -400.0), x, y);
            auto func = [](std::array<T, 2> arr) -> T {
                const auto& [x, y] = arr;
                return 400.0 * exp(sqr(x - 1.0) * -400.0);
            };
            CHECK_NOTATION("400.000000 x 1.000000 - 2.000000 ^ 400.000000 ~ * exp *");
            CHECK_VALUE(1., 100.);
            CHECK_VALUE(1.13, 95.);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR(exp(-t) * (sin(x)^2 + cos(x)^2), x, t);
            auto func = [](std::array<T, 2> arr) -> T {
                const auto& [x, t] = arr;
                return exp(-t) * (sqr(sin(x)) + sqr(cos(x)));
            };
            CHECK_NOTATION("t ~ exp x sin 2.000000 ^ x cos 2.000000 ^ + *");
            CHECK_VALUE(0.3, 1.7);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR(sqrt(x^2 + y^2) + abs(x - y) + (-sign(x*y)), x, y);
            auto func = [](std::array<T, 2> arr) -> T {
                const auto& [x, y] = arr;
                return sqrt(sqr(x) + sqr(y)) + abs(x - y) + (-sign(x*y));
            };
            CHECK_NOTATION("x 2.000000 ^ y 2.000000 ^ + sqrt x y - abs + x y * sign ~ +");
            CHECK_VALUE(3., 4.);
        }

        {
            SETUP_EXPR(atan2(y, x) + log1p(abs(x)) * expm1(y), x, y);
            CHECK_NOTATION("y x atan2 x abs log1p y expm1 * +");
            CHECK_VALUE(1., 0.5);
        }

        {
            SETUP_EXPR(abs(x) - y, x, y);
            CHECK_NOTATION("x abs y -");
            CHECK_VALUE(1., 5.5);
        }

        {
            SETUP_EXPR(atan2(y + 1, x * 2), x, y);
            CHECK_NOTATION("y 1.000000 + x 2.000000 * atan2");
            CHECK_VALUE(-1., 12.);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR((u + v - w)^3 / (1 + sqr(u - v)), u, v, w);
            auto func = [](std::array<T, 3> arr) -> T {
                const auto& [u, v, w] = arr;
                return pow(u + v - w, 3) / (1 + sqr(u - v));
            };
            CHECK_NOTATION("u v + w - 3.000000 ^ 1.000000 u v - sqr + /");
            CHECK_VALUE(2., 0.5, 1.);
        }

        {
            SETUP_EXPR(exp2(log2(x)) + expm1(log1p(x)), x);
            CHECK_NOTATION("x log2 exp2 x log1p expm1 +");
            CHECK_VALUE(0.25);
        }

        { // Operator '^' breaks the _FUNC macro, so we have to write this one manually
            PRESETUP_EXPR((a + b * c) / sqrt(a^2 + b^2 + c^2) - abs(a - b + c) * log(1 + exp(-a*b+c)), a, b, c);
            auto func = [](std::array<T, 3> arr) -> T {
                const auto& [a, b, c] = arr;
                return (a + b * c) / sqrt(sqr(a) + sqr(b) + sqr(c)) - abs(a - b + c) * log(1 + exp(-a*b+c));
            };
            CHECK_NOTATION("a b c * + a 2.000000 ^ b 2.000000 ^ + c 2.000000 ^ + sqrt / a b - c + abs 1.000000 a ~ b * c + exp + log * -");
            CHECK_VALUE(2., 0.5, 1.);
            CHECK_VALUE(1., 2, 3.);
            CHECK_VALUE(-1., 20, -3.);
        }

        {
            SETUP_EXPR(exp(sin(x) + cos(y)) * tan(x/y) - floor(abs(x) + ceil(y)) + sign(x*y) * atan2(x, y), x, y);
            CHECK_NOTATION("x sin y cos + exp x y / tan * x abs y ceil + floor - x y * sign x y atan2 * +");
            CHECK_VALUE(pi / 2., pi / 4.);
        }

        {
            PRESETUP_EXPR((p^q + r^s) / log(q + abs(r)) + min(p, max(q, r)) - sign(p*q*r*s) * (p - q + r - s)^3, p, q, r, s);
            auto func = [](std::array<T, 4> arr) -> T {
                const auto& [p, q, r, s] = arr;
                return (pow(p, q) + pow(r, s)) / log(q + abs(r)) + min(p, max(q, r)) - sign(p*q*r*s) * pow(p - q + r - s, 3);
            };
            CHECK_NOTATION("p q ^ r s ^ + q r abs + log / p q r max min + p q * r * s * sign p q - r + s - 3.000000 ^ * -");
            CHECK_VALUE(2., 3., 4., 5.);
        }

        {
            SETUP_EXPR(sqrt(sqr(u) * sqr(v)) + asin(u/v) - acos(v/u) * atan2(u, v) + (-sign(u + v)) * cbrt(u*v), u, v);
            CHECK_NOTATION("u sqr v sqr * sqrt u v / asin + v u / acos u v atan2 * - u v + sign ~ u v * cbrt * +");
            CHECK_VALUE(1., -1.);
        }

        {
            SETUP_EXPR(exp(m) * log(n) + sinh(m + n) - cosh(o) / tanh(m*n*o) + abs(fmod(m, n)) * round(o), m, n, o);
            CHECK_NOTATION("m exp n log * m n + sinh + o cosh m n * o * tanh / - m n fmod abs o round * +");
            CHECK_VALUE(0.5, 2., 3.);
        }

        {
            SETUP_EXPR((pow(a, b) - pow(b, a)) * sin(a*2/b) + cos(b*5/a) / tan(a + b) - sign(a - b) * sqrt(abs(a*b)), a, b);
            CHECK_NOTATION("a b pow b a pow - a 2.000000 * b / sin * b 5.000000 * a / cos a b + tan / + a b - sign a b * abs sqrt * -");
            CHECK_VALUE(4., 2.);
        }
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

#undef PRESETUP_EXPR
#undef SETUP_EXPR
#undef CHECK_VALUE
#undef CHECK_NOTATION
}
#include "utils.hpp"

#include <bitset>
#include <cctype>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>

namespace {
constexpr auto make_lut(std::string_view chars) {
    std::bitset<256> lut{};
    for (const auto c : chars) lut[static_cast<uint8_t>(c)] = true;
    return lut;
}
constexpr auto whitespace_chars = make_lut(" \t\n\r\v\f");
constexpr auto integer_chars = make_lut("0123456789");
constexpr auto float_chars = make_lut("0123456789.");
constexpr auto operator_chars = make_lut("+-*/^~");
constexpr auto symbol_chars = make_lut("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.");
constexpr auto separator_chars = make_lut(",:|");
constexpr auto parenthesis_chars = make_lut("()");
}  // namespace

namespace formula::utils {

std::ostream& operator<<(std::ostream& os, token_t token) {
    // ParenthesisRight - 16 chars
    os << '[' << std::setw(16);
    switch (token.type) {
        case token_t::type_t::Unknown:
            os << "Unknown";
            break;
        case token_t::type_t::Number:
            os << "Number";
            break;
        case token_t::type_t::ParenthesisLeft:
            os << "ParenthesisLeft";
            break;
        case token_t::type_t::ParenthesisRight:
            os << "ParenthesisRight";
            break;
        case token_t::type_t::Separator:
            os << "Separator";
            break;
        case token_t::type_t::Operator:
            os << "Operator";
            break;
        case token_t::type_t::Symbol:
            os << "Symbol";
            break;
    }
    os << "] : " << token.str;

    return os;
}

bool operator==(const token_t& lhs, const token_t& rhs) {
    return lhs.type == rhs.type && lhs.str == rhs.str;
}

std::vector<token_t> tokenize(std::string_view input) {
    constexpr std::string_view trim_chars = " \t\n\r\v\f";
    const auto first = input.find_first_not_of(trim_chars);
    if (first == std::string_view::npos) return {};
    const auto last = input.find_last_not_of(trim_chars);
    input = input.substr(first, last - first + 1);

    std::vector<token_t> tokens;

    enum class state_t : uint8_t {
        NewToken,
        IntegerNumber,
        RealNumber,
        ExponentNumber,
        Symbol,
        ParenthesisLeft,
        ParenthesisRight,
        Separator,
        Operator,
        CompleteToken,
        Terminate
    };
    state_t state = state_t::NewToken;

    // Scientific notation parsing context. Only meaningful in state_t::ExponentNumber.
    bool exponent_has_digit = false;
    bool exponent_sign_allowed = false;

    token_t::type_t token_type = token_t::type_t::Unknown;
    size_t parenthesis_checker = 0;
    std::string_view::const_iterator token_begin = input.end();

    for (auto it = input.begin(); state != state_t::Terminate;) {
        if(static_cast<uint8_t>(*it) > 127)
            throw std::domain_error{"Wrong expression format. The expression contains non-ANSI characters."};
        // If the iterator is at the end of the input we use '\0' as a sentinel value.
        const char c = (it == input.end()) ? '\0' : *it;
        // Handle incorrect implicit conversion. I think so... In any case, the app fails without this cast.
        uint8_t idx = static_cast<uint8_t>(c);

        switch (state) {
            case state_t::NewToken: {
                token_type = token_t::type_t::Unknown;
                token_begin = it;
                if (whitespace_chars[idx]) {
                    state = state_t::NewToken;
                } else if (float_chars[idx]) {
                    state = state_t::IntegerNumber;
                    continue;
                } else if (operator_chars[idx]) {
                    state = state_t::Operator;
                    continue;
                } else if (parenthesis_chars[idx]) {
                    state = c == '(' ? state_t::ParenthesisLeft : state_t::ParenthesisRight;
                    continue;
                } else if (separator_chars[idx]) {
                    state = state_t::Separator;
                    continue;
                } else {  // if nothing detected can be a symbol
                    state = state_t::Symbol;
                    continue;
                }
            } break;
            case state_t::IntegerNumber: {
                if (c == 'e' || c == 'E') {
                    exponent_has_digit = false;
                    exponent_sign_allowed = true;
                    state = state_t::ExponentNumber;
                } else if (!float_chars[idx] && symbol_chars[idx]) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid numeric construction"};
                } else if (!float_chars[idx]) {
                    token_type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    if (c == '.') state = state_t::RealNumber;
                }
            } break;
            case state_t::RealNumber: {
                if (c == 'e' || c == 'E') {
                    exponent_has_digit = false;
                    exponent_sign_allowed = true;
                    state = state_t::ExponentNumber;
                } else if (!integer_chars[idx] && symbol_chars[idx]) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid dots. Dots are only allowed in number representation"};
                } else if (!integer_chars[idx]) {
                    token_type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                }
            } break;
            case state_t::ExponentNumber: {
                // After 'e'/'E': optional sign, then require >= 1 digit.
                if (exponent_sign_allowed && (c == '+' || c == '-')) {
                    exponent_sign_allowed = false;
                } else if (integer_chars[idx]) {
                    exponent_has_digit = true;
                    exponent_sign_allowed = false;
                } else {
                    if (!exponent_has_digit)
                        throw std::domain_error{"Wrong expression format. The expression contains invalid exponent representation"};
                    if (symbol_chars[idx])
                        throw std::domain_error{"Wrong expression format. The expression contains invalid numeric construction"};

                    token_type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                }
            } break;
            case state_t::Operator: {
                token_type = token_t::type_t::Operator;
                state = state_t::CompleteToken;
            } break;
            case state_t::ParenthesisLeft: {
                ++parenthesis_checker;
                token_type = token_t::type_t::ParenthesisLeft;
                state = state_t::CompleteToken;
            } break;
            case state_t::ParenthesisRight: {
                --parenthesis_checker;
                token_type = token_t::type_t::ParenthesisRight;
                state = state_t::CompleteToken;
            } break;
            case state_t::Separator: {
                token_type = token_t::type_t::Separator;
                state = state_t::CompleteToken;
            } break;
            case state_t::Symbol: {
                if (!symbol_chars[idx]) {
                    token_type = token_t::type_t::Symbol;
                    state = state_t::CompleteToken;
                    continue;
                }
            } break;
            case state_t::CompleteToken: {
                if(token_begin == it)
                    throw std::runtime_error{"Wrong tokenizer behaviour. Empty token parsed."};
                tokens.emplace_back(token_t{token_type, {token_begin, it}});
                state = it == input.end() ? state_t::Terminate : state_t::NewToken;
                continue;
            } break;
        }
        if(it != input.end())
            ++it;
    }

    if (parenthesis_checker)
        throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
    return tokens;
}

std::unordered_map<std::string, std::size_t> get_variables(const std::span<token_t>& tokens) {
    std::unordered_map<std::string, std::size_t> variables;

    for (const auto& token : tokens) {
        if (token.type == token_t::type_t::Symbol) {
            variables[token.str] = variables.size();
            continue;
        }
        if (token.type == token_t::type_t::Separator && token.str == ",") {
            continue;
        }
        // numbers/operators/parentheses are invalid in variables declaration
        throw std::domain_error{"Wrong variables format. Variables must be declared as symbols before ':' separator."};
    }

    return variables;
}

const std::unordered_map<std::string, std::size_t>& get_operator_priority() {
    static const std::unordered_map<std::string, size_t> operator_priority{
        {"(", 0}, {"+", 1}, {"-", 1}, {"*", 2},
        {"/", 2}, {"^", 3}, {"pow", 3}, {"~", 4}, {"sin", 4}, 
        {"cos", 4}, {"tan", 4}, {"atan", 4}, {"exp", 4},
        {"abs", 4}, {"sign", 4}, {"sqr", 4},  {"sqrt", 4},
        {"log", 4}, {"tgamma", 4}, {"exp2", 4}, {"expm1", 4},
        {"log10", 4}, {"log2", 4}, {"log1p", 4}, {"cbrt", 4},
        {"asin", 4}, {"acos", 4}, {"sinh", 4}, {"cosh", 4},
        {"tanh", 4}, {"asinh", 4}, {"acosh", 4}, {"atanh", 4},
        {"erf", 4}, {"erfc", 4}, {"lgamma", 4}, {"ceil", 4},
        {"floor", 4}, {"round", 4}, {"trunc", 4}, {"atan2", 4}, 
        {"hypot", 4}, {"fmod", 4}, {"min", 4}, {"max", 4}
    };
    return operator_priority;
}


void check_variables_admissibility(const std::unordered_map<std::string, std::size_t>& variables) {
    const auto& operator_priority = get_operator_priority();
    for(const auto& [variable, _] : variables) 
        if (operator_priority.contains(variable))
            throw std::domain_error{"Invalid variable designation. Variable name <" + variable + "> \
                                        is unavailable."};   
}

}
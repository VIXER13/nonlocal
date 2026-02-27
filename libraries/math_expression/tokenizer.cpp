#include "tokenizer.hpp"

#include <array>
#include <cctype>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>

namespace {
constexpr auto make_lut(std::string_view chars) {
    std::array<bool, 256> lut{};
    for (const auto c : chars) lut[static_cast<uint8_t>(c)] = true;
    return lut;
}
constexpr auto whitespace_chars = make_lut(" \t\n\r\v\f");
constexpr auto integer_chars = make_lut("0123456789");
constexpr auto float_chars = make_lut("0123456789.");
constexpr auto operator_chars = make_lut("+-*/^");
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

std::vector<token_t> tokenize(std::string_view input) {
    input = input.substr(input.find_first_not_of(" "));
    input = input.substr(0, input.find_last_not_of(" ") + 1);
    if (input.empty()) return {};
    std::vector<token_t> tokens;

    enum class state_t : uint8_t {
        NewToken,
        IntegerNumber,
        RealNumber,
        Symbol,
        ParenthesisLeft,
        ParenthesisRight,
        Separator,
        Operator,
        CompleteToken,
        Terminate
    };
    state_t state = state_t::NewToken;

    token_t token;
    size_t parenthesis_checker = 0;

    for (auto it = input.begin(); state != state_t::Terminate;) {
        char c = (it == input.end() || *it > 255) ? 0 : *it;

        switch (state) {
            case state_t::NewToken: {
                // ToDo: Use iterators instead of string concatenation for better performance
                token = {token_t::type_t::Unknown, ""};
                if (whitespace_chars.at(c)) {
                    state = state_t::NewToken;
                } else if (float_chars.at(c)) {
                    state = state_t::IntegerNumber;
                    continue;
                } else if (operator_chars.at(c)) {
                    state = state_t::Operator;
                    continue;
                } else if (parenthesis_chars.at(c)) {
                    state = c == '(' ? state_t::ParenthesisLeft : state_t::ParenthesisRight;
                    continue;
                } else if (separator_chars.at(c)) {
                    state = state_t::Separator;
                    continue;
                } else {  // if nothing detected can be a symbol
                    state = state_t::Symbol;
                    continue;
                }
            } break;
            case state_t::IntegerNumber: {
                if (!float_chars.at(c) && symbol_chars.at(c)) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid numeric construction"};
                } else if (!float_chars.at(c)) {
                    token.type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    if (c == '.') state = state_t::RealNumber;
                    token.str.push_back(c);
                }
            } break;
            case state_t::RealNumber: {
                if (!integer_chars.at(c) && symbol_chars.at(c)) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid dots. Dots are only allowed in number representation"};
                } else if (!integer_chars.at(c)) {
                    token.type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    token.str.push_back(c);
                }
            } break;
            case state_t::Operator: {
                if (operator_chars.at(c)) {
                    token.str.push_back(c);
                    token.type = token_t::type_t::Operator;
                    state = state_t::CompleteToken;
                }
            } break;
            case state_t::ParenthesisLeft: {
                ++parenthesis_checker;
                token.str.push_back(c);
                token.type = token_t::type_t::ParenthesisLeft;
                state = state_t::CompleteToken;
            } break;
            case state_t::ParenthesisRight: {
                --parenthesis_checker;
                token.str.push_back(c);
                token.type = token_t::type_t::ParenthesisRight;
                state = state_t::CompleteToken;
            } break;
            case state_t::Separator: {
                token.str.push_back(c);
                token.type = token_t::type_t::Separator;
                state = state_t::CompleteToken;
            } break;
            case state_t::Symbol: {
                if (!symbol_chars.at(c)) {
                    token.type = token_t::type_t::Symbol;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    token.str.push_back(c);
                }
            } break;
            case state_t::CompleteToken: {
                tokens.push_back(std::move(token));
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
        // numbers/operators/parentheses are invalid in variables declaration
        throw std::domain_error{"Wrong variables format. Variables must be declared as symbols before ':' separator."};
    }

    return variables;
}
}  // namespace formula::utils

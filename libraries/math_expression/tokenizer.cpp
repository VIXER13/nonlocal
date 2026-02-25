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
constexpr auto operator_chars = make_lut("+-*/^");
constexpr auto symbol_chars = make_lut("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789.");
constexpr auto separator_chars = make_lut(",:|");
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

static std::vector<token_t> parse(std::string_view input) {
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
    };
    state_t state = state_t::NewToken;

    token_t token;
    size_t parenthesis_checker = 0;

    for (auto it = input.begin(); it != input.end() || state == state_t::CompleteToken;) {
        switch (state) {
            case state_t::NewToken: {
                token = {token_t::type_t::Unknown, ""};
                if (whitespace_chars.at(*it)) {
                    state = state_t::NewToken;
                } else if (integer_chars.at(*it)) {
                    state = state_t::IntegerNumber;
                    continue;
                } else if (operator_chars.at(*it)) {
                    state = state_t::Operator;
                    continue;
                } else if (*it == '(') {
                    state = state_t::ParenthesisLeft;
                    continue;
                } else if (*it == ')') {
                    state = state_t::ParenthesisRight;
                    continue;
                } else if (separator_chars.at(*it)) {
                    state = state_t::Separator;
                    continue;
                } else {  // if nothing detected can be a symbol
                    state = state_t::Symbol;
                    continue;
                }
            } break;
            case state_t::IntegerNumber: {
                if (!integer_chars.at(*it) && symbol_chars.at(*it)) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid numeric construction"};
                } else if (!integer_chars.at(*it)) {
                    token.type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    if (*it == '.') state = state_t::RealNumber;
                    token.str.push_back(*it);
                }
            } break;
            case state_t::RealNumber: {
                if (!integer_chars.at(*it) && symbol_chars.at(*it)) {
                    throw std::domain_error{"Wrong expression format. The expression contains invalid dots. Dots are only allowed in number representation"};
                } else if (!integer_chars.at(*it)) {
                    token.type = token_t::type_t::Number;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    token.str.push_back(*it);
                }
            } break;
            case state_t::Operator: {
                if (operator_chars.at(*it)) {
                    token.str.push_back(*it);
                    token.type = token_t::type_t::Operator;
                    state = state_t::CompleteToken;
                }
            } break;
            case state_t::ParenthesisLeft: {
                ++parenthesis_checker;
                token.str.push_back(*it);
                token.type = token_t::type_t::ParenthesisLeft;
                state = state_t::CompleteToken;
            } break;
            case state_t::ParenthesisRight: {
                --parenthesis_checker;
                token.str.push_back(*it);
                token.type = token_t::type_t::ParenthesisRight;
                state = state_t::CompleteToken;
            } break;
            case state_t::Separator: {
                token.str.push_back(*it);
                token.type = token_t::type_t::Separator;
                state = state_t::CompleteToken;
            } break;
            case state_t::Symbol: {
                if (!symbol_chars.at(*it)) {
                    token.type = token_t::type_t::Symbol;
                    state = state_t::CompleteToken;
                    continue;
                } else {
                    token.str.push_back(*it);
                }
            } break;
            case state_t::CompleteToken: {
                tokens.push_back(std::move(token));
                state = state_t::NewToken;
                continue;
            } break;
        }
        ++it;
    }

    if (parenthesis_checker)
        throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
    return tokens;
}

std::unordered_map<std::string, std::size_t> get_variables(const std::vector<token_t>& tokens) {
    std::unordered_map<std::string, std::size_t> variables;
    bool has_colon = false;

    for (const auto& token : tokens) {
        if (token.type == token_t::type_t::Separator) {
            if (token.str == ":") {
                has_colon = true;
                break;
            }
            continue;
        }

        else if (token.type == token_t::type_t::Symbol) {
            variables[token.str] = variables.size();
            continue;
        }

        // numbers/operators/parentheses before ':' are invalid in variables declaration
        throw std::domain_error{"Wrong variables format. Variables must be declared as symbols before ':' separator."};
    }

    if (!has_colon) {
        throw std::domain_error{"Wrong variables format. Symbol ':' is required after variables initialization."};
    }

    return variables;
}
}  // namespace formula::utils

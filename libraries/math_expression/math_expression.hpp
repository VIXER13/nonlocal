#pragma once

#include "utils.hpp"

#include <metamath/types/visitor.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <span>
#include <stack>
#include <stdexcept>
#include <variant>
#include <vector>

namespace formula {

template<class T>
class math_expression {
    struct variable_index final { size_t index; };
    using unary_operator = T(*)(const T);
    using binary_operator = T(*)(const T, const T);
    using operand = std::variant<T, variable_index, unary_operator, binary_operator>;

    std::vector<operand> _polish_notation;
    std::unordered_map<std::string, std::size_t> _variables;

#ifdef __clang__
    static inline std::stack<T, std::vector<T>> stack;
    static T pop_stack(std::stack<T, std::vector<T>>& stack);
#pragma omp threadprivate(stack)
#else
    static T pop_stack(std::stack<T>& stack);
#endif

    void assemble_polish_notation(const std::span<utils::token_t>& tokens);

    T calc_polish_notation(const std::span<const T> input_variables) const;

public:
    static const std::unordered_map<std::string, unary_operator>& unary_operators();
    static const std::unordered_map<std::string, binary_operator>& binary_operators();

    // pre_infix_notation -> |variables : expression|. Example: |x y z t : x * y * z - t / x + sin(x * y * z)|.
    // infix notation (standart) -> x + 5 * (y - z / t), polish notation (postfix) -> x 5 y z t / - * +
    explicit math_expression(std::string pre_infix_notation);

    std::string to_polish() const;
    std::size_t variables_count() const noexcept;
    const std::unordered_map<std::string, std::size_t>& variables() const noexcept;

    T operator()(const std::span<const T> input_vars) const;
    T operator()(const std::initializer_list<T>& input_vars) const;
};

template<class T>
const std::unordered_map<std::string, typename math_expression<T>::unary_operator>& math_expression<T>::unary_operators() {
    static constexpr auto unary_minus = [](const T value) noexcept { return -value; };
    static constexpr auto sign        = [](const T value) noexcept { return T((value > 0) - (value < 0)); };
    static constexpr auto sqr         = [](const T value) noexcept { return value * value; };
    static const std::unordered_map<std::string, typename math_expression<T>::unary_operator> operators{
        {"~", unary_minus}, {"sign", sign}, {"abs", std::abs}, {"sqr", sqr}, {"sqrt", std::sqrt}, {"cbrt", std::cbrt},
        {"sin", std::sin}, {"cos", std::cos}, {"tan", std::tan}, {"asin", std::asin}, {"acos", std::acos}, {"atan", std::atan},
        {"sinh", std::sinh}, {"cosh", std::cosh}, {"tanh", std::tanh}, {"asinh", std::asinh}, {"acosh", std::acosh}, {"atanh", std::atanh},
        {"exp", std::exp}, {"exp2", std::exp2}, {"expm1", std::expm1}, {"log", std::log}, {"log2", std::log2}, {"log1p", std::log1p}, {"log10", std::log10},
        {"erf", std::erf}, {"erfc", std::erfc}, {"tgamma", std::tgamma}, {"lgamma", std::lgamma},
        {"ceil", std::ceil}, {"trunc", std::trunc}, {"floor", std::floor}, {"round", std::round}
    };
    return operators;
}

template<class T>
const std::unordered_map<std::string, typename math_expression<T>::binary_operator>& math_expression<T>::binary_operators() {
    static constexpr auto plus       = [](const T left, const T right) noexcept { return left + right; };
    static constexpr auto minus      = [](const T left, const T right) noexcept { return left - right; };
    static constexpr auto divides    = [](const T left, const T right) noexcept { return left / right; };
    static constexpr auto multiplies = [](const T left, const T right) noexcept { return left * right; };
    static constexpr auto min        = [](const T left, const T right) noexcept { return std::min(left, right); };
    static constexpr auto max        = [](const T left, const T right) noexcept { return std::max(left, right); };
    static constexpr auto pow        = [](const T left, const T right) noexcept { return std::pow(left, right); };
    static const std::unordered_map<std::string, math_expression<T>::binary_operator> operators{
        {"+", plus}, {"-", minus}, {"*", multiplies}, {"/", divides}, {"^", std::pow}, {"pow", pow},
        {"atan2", std::atan2}, {"hypot", std::hypot}, {"fmod", std::fmod}, {"min", min}, {"max", max}
    };
    return operators;
}

template<class T>
void math_expression<T>::assemble_polish_notation(const std::span<utils::token_t>& tokens) {
    const auto& operator_priority = utils::get_operator_priority();
    const auto& unary = unary_operators();
    const auto& binary = binary_operators();
    std::vector<utils::token_t> prefixed;
    std::stack<utils::token_t> operators;
    auto prev_type = utils::token_t::type_t::Unknown;

    prefixed.reserve(tokens.size());

    for (auto& token: tokens) {
        switch(token.type) {
            case utils::token_t::type_t::Number:
            prefixed.push_back(token);
            prev_type = utils::token_t::type_t::Number;
            continue;
            case utils::token_t::type_t::ParenthesisLeft:
            operators.push(token);
            prev_type = utils::token_t::type_t::ParenthesisLeft;
            continue;
            case utils::token_t::type_t::ParenthesisRight: {
                while (!operators.empty() && operators.top().type != utils::token_t::type_t::ParenthesisLeft) {
                    prefixed.push_back(operators.top());
                    operators.pop();    
                }
                if (operators.empty())
                    throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
                operators.pop();

                // if a function (Symbol token) or unary operator is on top, pop it too
                if (!operators.empty() && operator_priority.contains(operators.top().str)
                    && (operators.top().type == utils::token_t::type_t::Symbol || operators.top().str == "~")) {
                    prefixed.push_back(operators.top());
                    operators.pop();
                }

                prev_type = utils::token_t::type_t::ParenthesisRight;
                continue;
            }
            case utils::token_t::type_t::Symbol:
            case utils::token_t::type_t::Operator: {
                if (operator_priority.contains(token.str)) {
                    if(token.str == "-") {
                        switch(prev_type) {
                            case utils::token_t::type_t::Unknown:
                            case utils::token_t::type_t::Operator:
                            case utils::token_t::type_t::Separator:
                            case utils::token_t::type_t::ParenthesisLeft:
                                token.str = "~";
                            default:
                                break;
                        }
                    }

                    const bool token_is_right_associative = (token.str == "^") || token.str == "pow" || unary.contains(token.str);

                    while (!operators.empty() && operators.top().type != utils::token_t::type_t::ParenthesisLeft) {
                        const auto top_priority = operator_priority.at(operators.top().str);
                        const auto token_priority = operator_priority.at(token.str);
                        const bool should_pop = (top_priority > token_priority) || (top_priority == token_priority && !token_is_right_associative);

                        if (!should_pop)
                            break;

                        prefixed.push_back(operators.top());
                        operators.pop();
                    }
                    operators.push(token);
                    prev_type = utils::token_t::type_t::Operator;
                    continue;
                } else if (_variables.contains(token.str)) {
                    prefixed.push_back(token);
                    prev_type = utils::token_t::type_t::Symbol;
                    continue;
                }

                throw std::domain_error{"Wrong variables and operators searching. Can not find contained string."};
            }
            case utils::token_t::type_t::Separator: {
                if (token.str != ",")
                    throw std::domain_error{"Wrong expression format. Unexpected separator."};
                while (!operators.empty() && operators.top().type != utils::token_t::type_t::ParenthesisLeft) {
                    prefixed.push_back(operators.top());
                    operators.pop();
                }
                if (operators.empty())
                    throw std::domain_error{"Wrong expression format. Separator used out of parentheses."};
                prev_type = utils::token_t::type_t::Separator;
                continue;
            }
            default: {
                throw std::domain_error{"Wrong expression format. Unexpected token."};
            }
        }
        throw std::domain_error{"Wrong expression format. Unexpected token."};
    }

    while (!operators.empty()) {
        if (operators.top().type == utils::token_t::type_t::ParenthesisLeft)
            throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
        prefixed.push_back(operators.top());
        operators.pop();
    }

    _polish_notation.reserve(prefixed.size());
    for(const auto& token : prefixed) {
        switch(token.type) {
            case utils::token_t::type_t::Number:
                _polish_notation.push_back(utils::get_number<T>(token.str));
                continue;
            case utils::token_t::type_t::Symbol:
            case utils::token_t::type_t::Operator:
                if (token.type == utils::token_t::type_t::Symbol && _variables.contains(token.str))
                    _polish_notation.push_back(variable_index{_variables.at(token.str)});
                else if (binary.contains(token.str))
                    _polish_notation.push_back(binary.at(token.str));
                else if(unary.contains(token.str))
                    _polish_notation.push_back(unary.at(token.str));
                else
                    throw std::domain_error{"Wrong variables and operators searching. Can not find contained string."};
                continue;
            default:
                throw std::domain_error{"Wrong expression format. Unexpected token."};
        }
    }
}

template<class T>
#ifdef __clang__
T math_expression<T>::pop_stack(std::stack<T, std::vector<T>>& stack) {
#else
T math_expression<T>::pop_stack(std::stack<T>& stack) {
#endif
    const T value = stack.top();
    stack.pop();
    return value;
}

template<class T>
T math_expression<T>::calc_polish_notation(const std::span<const T> input_variables) const {
    if (input_variables.size() != _variables.size())
        throw std::domain_error{"Wrong number of variables."};
#ifndef __clang__
    std::stack<T> stack;
#endif
    for(const auto& operation : _polish_notation)
        std::visit(metamath::types::visitor{
            [&](const T number) { stack.push(number); },
            [&](const variable_index& index) { stack.push(input_variables[index.index]); },
            [&](const unary_operator& operation) { stack.push(operation(pop_stack(stack))); },
            [&](const binary_operator& operation) {
                const T right = pop_stack(stack);
                const T left = pop_stack(stack);
                stack.push(operation(left, right));
            }
        }, operation);
    return pop_stack(stack);
}

template<class T>
math_expression<T>::math_expression(std::string pre_infix_notation) {
    auto tokens = utils::tokenize(pre_infix_notation);
    auto sep_it = std::find_if(tokens.begin(), tokens.end(), [](const utils::token_t& t) {
        return t.type == utils::token_t::type_t::Separator && t.str == ":";
    });
    if (sep_it == tokens.end())
        throw std::domain_error{"Wrong variables format. Symbol ':' is required after variables initialization."};
    _variables = utils::get_variables(std::span(tokens.begin(), sep_it));
    utils::check_variables_admissibility(_variables);
    auto infix_notation = std::span(std::next(sep_it), tokens.end());
    if (infix_notation.empty())
        throw std::domain_error{"Wrong expression format. Formula after ':' is required."};
    assemble_polish_notation(infix_notation);
}

template<class T>
std::string math_expression<T>::to_polish() const {
    static constexpr auto get_name = [](const auto& operation, const auto& values) {
        for (const auto& [name, index] : values)
            if (index == operation)
                return name;
        throw std::domain_error{"Unknown operator, unable to recover polish notation."};
    };
    const auto concatenate = [this](const std::string& polish, const operand& operation) {
        return polish + (polish.empty() ? "" : " ") + std::visit(metamath::types::visitor{
            [](const T number) { return std::to_string(number); },
            [this](const variable_index& index) { return get_name(index.index, variables()); },
            [this](const unary_operator& operation) { return get_name(operation, unary_operators()); },
            [this](const binary_operator& operation) { return get_name(operation, binary_operators()); },
        }, operation);
    };
    using namespace std::string_literals;
    return std::accumulate(_polish_notation.begin(), _polish_notation.end(), ""s, concatenate);
}

template<class T>
std::size_t math_expression<T>::variables_count() const noexcept {
    return _variables.size();
}

template<class T>
const std::unordered_map<std::string,std::size_t>& math_expression<T>::variables() const noexcept {
    return _variables;
}

template<class T>
T math_expression<T>::operator()(const std::span<const T> input_vars) const {
    return calc_polish_notation(input_vars);
}

template<class T>
T math_expression<T>::operator()(const std::initializer_list<T>& input_vars) const {
    return (*this)(std::span(input_vars));
}

}
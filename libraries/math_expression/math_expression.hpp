#pragma once

#include "tokenizer.hpp"
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

    static const std::unordered_map<std::string, unary_operator>& unary_operators();
    static const std::unordered_map<std::string, binary_operator>& binary_operators();

    std::unordered_map<std::size_t, std::string> find_variables_and_operators(const std::string& infix_notation) const;
    void assemble_polish_notation(const std::string& infix_notation);
    void assemble_polish_notation(const std::vector<utils::token_t>& tokens);

    T calc_polish_notation(const std::span<const T> input_variables) const;

public:
    // pre_infix_notation -> |variables : expression|. Example: |x y z t : x * y * z - t / x + sin(x * y * z)|.
    // infix notation (standart) -> x + 5 * (y - z / t), polish notation (prefix) -> x 5 y z t / - * +
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
    static const std::unordered_map<std::string, math_expression<T>::binary_operator> operators{
        {"+", plus}, {"-", minus}, {"*", multiplies}, {"/", divides}, {"^", std::pow}
    };
    return operators;
}

template<class T>
std::unordered_map<std::size_t, std::string> math_expression<T>::find_variables_and_operators(const std::string& infix_notation) const {
    const auto& operator_priority = utils::get_operator_priority();
    const auto& one_symbol_operator = utils::get_one_sym_operators();
    std::unordered_map<std::size_t, std::string> variables_and_operators_indices;

    const auto check_integrity = [&one_symbol_operator, &infix_notation](std::size_t pos, const std::string& smth) {
        const std::size_t shift = smth.size() - 1;
        const bool is_not_boundary_position = ((0 < pos && pos < infix_notation.size() - 1 - shift)
            && (one_symbol_operator.contains(infix_notation[pos - 1]) && one_symbol_operator.contains(infix_notation[pos + shift + 1]) ));
        const bool is_boundary_position = (pos == 0 && one_symbol_operator.contains(infix_notation[pos + shift + 1]) ) ||
            (pos == infix_notation.size() - 1 - shift && one_symbol_operator.contains(infix_notation[pos - 1]) );
        return is_boundary_position || is_not_boundary_position;
    };

    const auto push_it = [&infix_notation, &check_integrity, &variables_and_operators_indices](const std::string& smth) {
        for(std::size_t position = infix_notation.find(smth); position != std::string::npos; position = infix_notation.find(smth, position + smth.size()))
            if (check_integrity(position, smth))
                variables_and_operators_indices[position] = smth;
    };

    for(const auto& [variable, _] : _variables)
        push_it(variable);

    for(const auto& [op, priority] : operator_priority)
        if (!one_symbol_operator.contains(op.front()))
            push_it(op);

    return variables_and_operators_indices;
}

template<class T>
void math_expression<T>::assemble_polish_notation(const std::string& infix_notation) {
    const auto variables_and_operators_indices = find_variables_and_operators(infix_notation);
    const auto& operator_priority = utils::get_operator_priority();
    const auto& one_symbol_operator = utils::get_one_sym_operators();
    std::stack<std::string> operators;
    std::vector<std::string> polish_notation;
    for(std::size_t i = 0; i < infix_notation.size(); ++i) {
        const char symbol = infix_notation[i];
        if (variables_and_operators_indices.contains(i)) {
            const auto& smth = variables_and_operators_indices.at(i);
            if (operator_priority.contains(smth)) {
                while (!operators.empty() && (operator_priority.at(operators.top()) >= operator_priority.at(smth))) {
                    polish_notation.push_back(operators.top());
                    operators.pop();
                }
                operators.push(smth);
            } else if (_variables.contains(smth)) {
                polish_notation.push_back(smth);
            } else {
                throw std::domain_error{"Wrong variables and operators searching. Can not find contained string."};
            }
            i += smth.size() - 1;
        } else if (std::isdigit(symbol) || symbol == '.') {
            if (i == 0 || (i > 0 && !std::isdigit(infix_notation[i - 1]) && infix_notation[i - 1] != '.')) {
                polish_notation.emplace_back(std::string{ symbol });
            } else {
                polish_notation.back().push_back(symbol);
            }
        } else if (symbol == '(') {
            operators.push(std::string{ symbol });
        } else if (symbol == ')') {
            while (!operators.empty() && operators.top() != std::string{ '(' }) {
                polish_notation.push_back(operators.top());
                operators.pop();
            }
            operators.pop();
        } else if (one_symbol_operator.contains(symbol) && operator_priority.contains(std::string{symbol})) {
            std::string op = std::string{ symbol };
            // for unary minus
            if (op == std::string{ '-' } && (i == 0 || (i > 1 && operator_priority.contains(std::string{ infix_notation[i - 1] }) )))
                op = std::string{ '~' };

            while (!operators.empty() && (operator_priority.at(operators.top()) >= operator_priority.at(op))) {
                polish_notation.push_back(operators.top());
                operators.pop();
            }
            operators.push(op);
        }
    }
    while (!operators.empty()) {
        polish_notation.push_back(operators.top());
        operators.pop();
    }

    const auto& unary = unary_operators();
    const auto& binary = binary_operators();
    _polish_notation.reserve(polish_notation.size());
    for(const std::string& smth : polish_notation) {
        if (utils::is_number(smth))
            _polish_notation.push_back(utils::get_number<T>(smth));
        else if (_variables.count(smth))
            _polish_notation.push_back(variable_index{_variables.at(smth)});
        else if (binary.count(smth))
            _polish_notation.push_back(binary.at(smth));
        else
            _polish_notation.push_back(unary.at(smth));
    }
}

template<class T>
void math_expression<T>::assemble_polish_notation(const std::vector<utils::token_t>& tokens) {
    const auto& operator_priority = utils::get_operator_priority();

    // expression starts after ':'
    std::size_t expr_start = 0;
    for (; expr_start < tokens.size(); ++expr_start) {
        if (tokens[expr_start].type == utils::token_t::type_t::Separator && tokens[expr_start].str == ":") {
            ++expr_start;
            break;
        }
    }
    if (expr_start == 0 || expr_start > tokens.size())
        throw std::domain_error{"Wrong variables format. Symbol ':' is required after variables initialization."};

    const auto& unary = unary_operators();
    const auto& binary = binary_operators();

    using op_entry = std::pair<int, operand>; // (priority, operator operand). priority < 0 means '('
    std::stack<op_entry> operators;

    const auto is_paren = [](const op_entry& e) { return e.first < 0; };

    const auto emit_operator = [this, &is_paren](const op_entry& e) {
        if (is_paren(e))
            throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
        _polish_notation.push_back(e.second);
    };

    const auto make_paren = []() -> op_entry {
        // sentinel: '(' marker
        return {-1, operand{unary_operator{nullptr}}};
    };

    const auto make_operator = [&operator_priority, &unary, &binary](const std::string& op) -> op_entry {
        if (!operator_priority.contains(op))
            throw std::domain_error{"Wrong expression format. Unknown operator."};
        const int priority = static_cast<int>(operator_priority.at(op));
        if (binary.contains(op))
            return {priority, operand{binary.at(op)}};
        if (unary.contains(op))
            return {priority, operand{unary.at(op)}};
        throw std::domain_error{"Wrong expression format. Unknown operator."};
    };

    const auto is_value_token = [this](const utils::token_t& t) {
        return t.type == utils::token_t::type_t::Number ||
               (t.type == utils::token_t::type_t::Symbol && _variables.contains(t.str)) ||
               t.type == utils::token_t::type_t::ParenthesisRight;
    };

    _polish_notation.clear();
    _polish_notation.reserve(tokens.size() - expr_start);

    for (std::size_t i = expr_start; i < tokens.size(); ++i) {
        const auto& token = tokens[i];

        if (token.type == utils::token_t::type_t::Separator) {
            // allow optional wrappers like '|'
            if (token.str == "|" || token.str == ",")
                continue;
            throw std::domain_error{"Wrong expression format. Unexpected separator in expression."};
        }

        if (token.type == utils::token_t::type_t::Number) {
            _polish_notation.push_back(utils::get_number<T>(token.str));
            continue;
        }

        if (token.type == utils::token_t::type_t::ParenthesisLeft) {
            operators.push(make_paren());
            continue;
        }

        if (token.type == utils::token_t::type_t::ParenthesisRight) {
            while (!operators.empty() && !is_paren(operators.top())) {
                emit_operator(operators.top());
                operators.pop();
            }
            if (operators.empty())
                throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
            operators.pop();
            continue;
        }

        if (token.type == utils::token_t::type_t::Symbol) {
            if (_variables.contains(token.str)) {
                _polish_notation.push_back(variable_index{_variables.at(token.str)});
                continue;
            }

            // function-like operators: sin, cos, ... are tokenized as Symbol
            if (operator_priority.contains(token.str)) {
                const op_entry current = make_operator(token.str);
                while (!operators.empty() && !is_paren(operators.top()) && operators.top().first >= current.first) {
                    emit_operator(operators.top());
                    operators.pop();
                }
                operators.push(current);
                continue;
            }

            throw std::domain_error{"Wrong variables and operators searching. Can not find contained string."};
        }

        if (token.type == utils::token_t::type_t::Operator) {
            std::string op = token.str;

            // unary minus
            if (op == "-") {
                const bool unary_context = (i == expr_start) || (i > expr_start && !is_value_token(tokens[i - 1]));
                if (unary_context)
                    op = "~";
            }

            const op_entry current = make_operator(op);
            while (!operators.empty() && !is_paren(operators.top()) && operators.top().first >= current.first) {
                emit_operator(operators.top());
                operators.pop();
            }
            operators.push(current);
            continue;
        }

        throw std::domain_error{"Wrong expression format. Unexpected token."};
    }

    while (!operators.empty()) {
        if (is_paren(operators.top()))
            throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
        emit_operator(operators.top());
        operators.pop();
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
    const std::size_t delimiter = pre_infix_notation.find(':');
    if (delimiter == std::string::npos)
        throw std::domain_error{"Wrong variables format. Symbol ':' is required after variables initialization."};
    _variables = utils::get_variables(utils::trim(pre_infix_notation.substr(0, delimiter)));
    utils::check_variables_admissibility(_variables);
    std::string infix_notation = pre_infix_notation.substr(delimiter + 1);
    infix_notation = utils::delete_all(infix_notation, ' ');
    if (infix_notation.empty())
        throw std::domain_error{"Wrong expression format. Formula after ':' is required."};
    utils::parentheses_check(infix_notation);
    utils::dots_check(infix_notation);
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
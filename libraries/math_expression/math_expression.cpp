#include "math_expression.hpp"

#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <iostream>
#include <ranges>

namespace {
const std::unordered_map<std::string, std::size_t>& get_operator_priority() {
    using namespace std::string_literals;
    static const std::unordered_map<std::string, std::size_t> operator_priority{{"("s, 0}, {"+"s, 1}, {"-"s, 1}, {"*"s, 2},
                                                                                {"/"s, 2}, {"^"s, 3}, {"~"s, 4}, {"sin"s, 4}, 
                                                                                {"cos"s, 4}, {"tan"s, 4}, {"atan"s, 4}, {"exp"s, 4},
                                                                                {"abs"s, 4}, {"sign"s, 4}, {"sqr"s, 4},  {"sqrt"s, 4},
                                                                                {"log"s, 4}, {"tgamma"s, 4}, {"exp2"s, 4}, {"expm1"s, 4},
                                                                                {"log10"s, 4}, {"log2"s, 4}, {"log1p"s, 4}, {"cbrt"s, 4},
                                                                                {"asin"s, 4}, {"acos"s, 4}, {"sinh"s, 4}, {"cosh"s, 4},
                                                                                {"tanh"s, 4}, {"asinh"s, 4}, {"acosh"s, 4}, {"atanh"s, 4},
                                                                                {"erf"s, 4}, {"erfc"s, 4}, {"lgamma"s, 4}, {"ceil"s, 4},
                                                                                {"floor"s, 4}, {"round"s, 4}, {"trunc"s, 4}};
    return operator_priority;
}

const std::unordered_set<char>& get_one_sym_operators() {
    static const std::unordered_set<char> one_sym_operators {'(', ')', '+', '-', '*', '/', '~', '^'};
    return one_sym_operators;
}

void check_variables_admissibility(const std::unordered_map<std::string, std::size_t>& variables) {
    const auto& operator_priority = get_operator_priority();
    for (const auto& [variable, _] : variables) 
        if (operator_priority.contains(variable))
            throw std::domain_error{"Invalid variable designation. Variable name <" + variable + "> \
                                        is unavailable."};   
}

void parentheses_check(const std::string& infix_notation) {
    if (formula::utils::count_all(infix_notation, '(') != formula::utils::count_all(infix_notation, ')'))
        throw std::domain_error{"Wrong expression format. The expression contains open parentheses."};
}

void dots_check(const std::string& infix_notation) {
    static constexpr char dot = '.';
    if (infix_notation[0] == dot || infix_notation[infix_notation.size() - 1] == dot)
        throw std::domain_error{"Wrong expression format. The expression can not be started or finished with dot. Dots are only allowed in number representation."};
    for (std::size_t position = infix_notation.find(dot); position != std::string::npos; position = infix_notation.find(dot, position + 1)) {
        if (!std::isdigit(infix_notation[position + 1]))
            throw std::domain_error{"Wrong expression format. The expression contains invalid dots. Dots are only allowed in number representation"};
    }
}

void variables_format_check(const std::vector<std::string>& variables) {
    for (const std::string& var : variables) 
        if (!var.empty() && !std::isalpha(var.front()))
            throw std::domain_error{"Wrong variables format <" + var + ">.\
                                        Variables must start with latin letter, variable cannot start with a number."};
}

std::unordered_map<std::string, std::size_t> get_variables(std::string pre_variables) {
    static constexpr char delimiter_symbol = ' ';
    std::vector<std::string> variables = {};
    std::size_t delimiter = pre_variables.find(delimiter_symbol);   
    while (delimiter != std::string::npos || pre_variables.size() > 0) {
        variables.push_back(pre_variables.substr(0, delimiter));
        pre_variables = delimiter < pre_variables.size() ? pre_variables.substr(delimiter + 1) : "";
        delimiter = pre_variables.find(delimiter_symbol);
    }
    variables_format_check(variables);
    std::unordered_map<std::string, std::size_t> res;
    std::ranges::for_each(std::views::iota(size_t(0), variables.size()), [&res, &variables](std::size_t i) { 
            res[variables[i]] = i; 
    });
    return res;
}

}

namespace formula {

math_expression::math_expression(std::string pre_infix_notation) {
    const std::size_t delimiter = pre_infix_notation.find(':');
    if (delimiter == std::string::npos)
        throw std::domain_error{"Wrong variables format. Symbol ':' is required after variables initialization."};
    _variables = ::get_variables(utils::trim(pre_infix_notation.substr(0, delimiter)));
    check_variables_admissibility(_variables);
    std::string infix_notation = pre_infix_notation.substr(delimiter + 1);
    infix_notation = utils::delete_all(infix_notation, ' ');
    if (infix_notation.empty())
        throw std::domain_error{"Wrong expression format. Formula after ':' is required."};
    parentheses_check(infix_notation);
    dots_check(infix_notation);
    assemble_polish_notation(infix_notation);
}

std::string math_expression::to_polish() const {
    auto concatenate = [](const std::string& a, const std::string& b){ return a + b; };
    return std::accumulate(_polish_notation.begin(), _polish_notation.end(), std::string{""}, concatenate);
}

std::size_t math_expression::variables_count() const {
    return _variables.size();
}

std::unordered_map<std::string,std::size_t> math_expression::get_variables() const {
    return _variables;
}

const std::unordered_map<std::string, math_expression::operator_index>& math_expression::get_arithmetic_operators() {
    using namespace std::string_literals;
    static const std::unordered_map<std::string, operator_index>
        arithmetic_operators{  {"+"s,     operator_index::plus},   {"-"s,      operator_index::minus},  {"*"s,      operator_index::multiply},
                               {"/"s,     operator_index::divide}, {"^"s,      operator_index::power},  {"~"s,      operator_index::unary_minus}, 
                               {"sin"s,   operator_index::sin},    {"cos"s,    operator_index::cos},    {"tan"s,    operator_index::tan},  
                               {"atan"s,  operator_index::atan},   {"exp"s,    operator_index::exp},    {"abs"s,    operator_index::abs},   
                               {"sign"s,  operator_index::sign},   {"sqr"s,    operator_index::sqr},    {"sqrt"s,   operator_index::sqrt},
                               {"log"s,   operator_index::log},    {"tgamma"s, operator_index::tgamma}, {"exp2"s,   operator_index::exp2},
                               {"expm1"s, operator_index::expm1},  {"log10"s,  operator_index::log10},  {"log2"s,   operator_index::log2},
                               {"log1p"s, operator_index::log1p},  {"cbrt"s,   operator_index::cbrt},   {"asin"s,   operator_index::asin},
                               {"acos"s,  operator_index::acos},   {"sinh"s,   operator_index::sinh},   {"cosh"s,   operator_index::cosh},
                               {"tanh"s,  operator_index::tanh},   {"asinh"s,  operator_index::asinh},  {"acosh"s,  operator_index::acosh},
                               {"atanh"s, operator_index::atanh},  {"erf"s,    operator_index::erf},    {"lgamma"s, operator_index::lgamma},
                               {"erfc"s,  operator_index::erfc},   {"ceil"s,   operator_index::ceil},   {"trunc"s,  operator_index::trunc},
                               {"floor"s, operator_index::floor},  {"round"s,  operator_index::round}};
    return arithmetic_operators;
}

std::unordered_map<std::size_t, std::string> math_expression::find_variables_and_operators(const std::string& infix_notation) const {
    const auto& operator_priority = get_operator_priority();
    const auto& one_symbol_operator = get_one_sym_operators();
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
        for (std::size_t position = infix_notation.find(smth); position != std::string::npos; position = infix_notation.find(smth, position + smth.size()))
            if (check_integrity(position, smth))
                variables_and_operators_indices[position] = smth;
    };

    for (const auto& [variable, _] : _variables)
        push_it(variable);

    for (const auto& [op, priority] : operator_priority)
        if (!one_symbol_operator.contains(op.front()))
            push_it(op);

    if (variables_and_operators_indices.empty()) 
        std::cout << "Warning: expression does not depend on the variables." << std::endl;
    return variables_and_operators_indices;
}

void math_expression::assemble_polish_notation(const std::string& infix_notation) {
    const auto variables_and_operators_indices = find_variables_and_operators(infix_notation);
    const auto& operator_priority = get_operator_priority();
    const auto& one_symbol_operator = get_one_sym_operators();
    std::stack<std::string> operators;
    for (std::size_t i = 0; i < infix_notation.size(); ++i) {
        const char symbol = infix_notation[i];
        if (variables_and_operators_indices.contains(i)) {
            const auto& smth = variables_and_operators_indices.at(i);
            if (operator_priority.contains(smth)) {
                while (!operators.empty() && (operator_priority.at(operators.top()) >= operator_priority.at(smth))) {
                    _polish_notation.push_back(operators.top());
                    operators.pop();
                }
                operators.push(smth);
            } else if (_variables.contains(smth)) {
                _polish_notation.push_back(smth);
            } else {
                throw std::domain_error{"Wrong variables and operators searching. Can not find contained string."};
            }
            i += smth.size() - 1;
        } else if (std::isdigit(symbol) || symbol == '.') {
            if (i == 0 || (i > 0 && !std::isdigit(infix_notation[i - 1]) && infix_notation[i - 1] != '.')) {
                _polish_notation.emplace_back(std::string{ symbol });
            } else {
                _polish_notation.back().push_back(symbol);
            }
        } else if (symbol == '(') {
            operators.push(std::string{ symbol });
        } else if (symbol == ')') {
            while (!operators.empty() && operators.top() != std::string{ '(' }) {
                _polish_notation.push_back(operators.top());
                operators.pop();
            }
            operators.pop();
        } else if (one_symbol_operator.contains(symbol) && operator_priority.contains(std::string{symbol})) {
            std::string op = std::string{ symbol };
            // for unary minus
            if (op == std::string{ '-' } && (i == 0 || (i > 1 && operator_priority.contains(std::string{ infix_notation[i - 1] }) )))
                op = std::string{ '~' };

            while (!operators.empty() && (operator_priority.at(operators.top()) >= operator_priority.at(op))) {
                _polish_notation.push_back(operators.top());
                operators.pop();
            }
            operators.push(op);
        }
    }
    while (!operators.empty()) {
        _polish_notation.push_back(operators.top());
        operators.pop();
    }
}

}; 
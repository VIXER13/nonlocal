#include "utils.hpp"

#include <algorithm>
#include <numeric>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace {

void variables_format_check(const std::vector<std::string>& variables) {
    for(const std::string& var : variables) 
        if (!var.empty() && !std::isalpha(var.front()))
            throw std::domain_error{"Wrong variables format <" + var + ">. "
                                    "Variables must start with latin letter, variable cannot start with a number."};
}

}

namespace formula::utils {

std::string trim(std::string str)
{
    str.erase(str.find_last_not_of(' ')+1);    //suffixing spaces
    str.erase(0, str.find_first_not_of(' '));  //prefixing spaces
    return str;
}

std::string& delete_all(std::string& str, char symb) {
    str.erase(std::remove(str.begin(), str.end(), symb), str.end());
    return str;
}

std::size_t count_all(const std::string& s, char symb) {
    return std::accumulate(s.begin(), s.end(), std::size_t(0), [&](std::size_t res, char a){ return res + std::size_t(a == symb); });
}

bool is_latin_str(const std::string& s) {
    if (s.size() == 0 || s.size() > 1)
        throw std::domain_error{"Wrong format, latyn symbol contains of one symbol."};
    for(const char d : s)
        if (!std::isalpha(d)) 
            return false;
    return true;
};

bool is_number(const std::string& s) {
    if (const std::size_t dots = count_all(s, '.'); dots > 1)
        return false;
    for(const char d : s)
        if (d != '.' && !std::isdigit(d)) 
            return false;
    return true;
};

const std::unordered_map<std::string, std::size_t>& get_operator_priority() {
    static const std::unordered_map<std::string, size_t> operator_priority{
        {"(", 0}, {"+", 1}, {"-", 1}, {"*", 2},
        {"/", 2}, {"^", 3}, {"~", 4}, {"sin", 4}, 
        {"cos", 4}, {"tan", 4}, {"atan", 4}, {"exp", 4},
        {"abs", 4}, {"sign", 4}, {"sqr", 4},  {"sqrt", 4},
        {"log", 4}, {"tgamma", 4}, {"exp2", 4}, {"expm1", 4},
        {"log10", 4}, {"log2", 4}, {"log1p", 4}, {"cbrt", 4},
        {"asin", 4}, {"acos", 4}, {"sinh", 4}, {"cosh", 4},
        {"tanh", 4}, {"asinh", 4}, {"acosh", 4}, {"atanh", 4},
        {"erf", 4}, {"erfc", 4}, {"lgamma", 4}, {"ceil", 4},
        {"floor", 4}, {"round", 4}, {"trunc", 4}
    };
    return operator_priority;
}

const std::unordered_set<char>& get_one_sym_operators() {
    static const std::unordered_set<char> one_sym_operators{'(', ')', '+', '-', '*', '/', '~', '^'};
    return one_sym_operators;
}

void check_variables_admissibility(const std::unordered_map<std::string, std::size_t>& variables) {
    const auto& operator_priority = get_operator_priority();
    for(const auto& [variable, _] : variables) 
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
    for(std::size_t position = infix_notation.find(dot); position != std::string::npos; position = infix_notation.find(dot, position + 1)) {
        if (!std::isdigit(infix_notation[position + 1]))
            throw std::domain_error{"Wrong expression format. The expression contains invalid dots. Dots are only allowed in number representation"};
    }
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
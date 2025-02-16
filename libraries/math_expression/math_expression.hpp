#pragma once

#include "utils.hpp"

#include <unordered_map>
#include <stack>
#include <vector>
#include <span>

#include <cmath>
#include <stdexcept>

namespace formula {

class math_expression {
private:
    std::vector<std::string> _polish_notation{};
    std::unordered_map<std::string, std::size_t> _variables;

public:
    // pre_infix_notation -> |variables : expression|. Example: |x y z t : x * y * z - t / x + sin(x * y * z)|.
    // infix notation (standart) -> x + 5 * (y - z / t), polish notation (prefix) -> x 5 y z t / - * +
    explicit math_expression(std::string pre_infix_notation);

    std::string to_polish() const;
    std::size_t variables_count() const;
    std::unordered_map<std::string, std::size_t> get_variables() const;

    template <utils::arithmetic T>
    T operator()(const std::span<const T> input_vars) const {
        return calc_polish_notation(input_vars);
    }

    template <utils::arithmetic T>
    T operator()(const std::initializer_list<T>& input_vars) const {
        return this->operator()(std::span(input_vars));
    }

private:
    enum class operator_index : std::size_t {
        plus, minus, unary_minus,
        multiply, divide, power, sqr, sqrt, cbrt,
        sin, asin, sinh, asinh, cos, acos, cosh, acosh, tan, atan, tanh, atanh,
        exp, exp2, expm1, log, log10, log2, log1p,
        abs, sign, ceil, floor, trunc, round,
        tgamma, lgamma, erf, erfc
    };

    template<utils::arithmetic T>
    T execute(const std::unordered_map<std::string, math_expression::operator_index>& ops, 
              std::stack<T>& calculation_values, const std::string& op) const;

    template <utils::arithmetic T>
    T calc_polish_notation(const std::span<const T> input_variables) const;

    static const std::unordered_map<std::string, operator_index>& get_arithmetic_operators();
    std::unordered_map<std::size_t, std::string> find_variables_and_operators(const std::string& infix_notation) const;
    void assemble_polish_notation(const std::string& infix_notation);
};

template<utils::arithmetic T>
T math_expression::execute(const std::unordered_map<std::string, math_expression::operator_index>& ops, 
          std::stack<T>& calculation_values, const std::string& op) const {
    const auto pop_element = [&calculation_values]() {
        T elem = T(0);
        if (!calculation_values.empty()) {
            elem = calculation_values.top();
            calculation_values.pop();
        }
        return elem;
    };
    const T right = pop_element();
    switch(ops.at(op))
    {
    case operator_index::plus:
        return pop_element() + right;
    case operator_index::minus:
        return pop_element() - right;
    case operator_index::multiply:
        return pop_element() * right;
    case operator_index::divide:
        return pop_element() / right;
    case operator_index::power:
        return std::pow(pop_element(), right);
    case operator_index::unary_minus:
        return -right;
    case operator_index::sin:
        return std::sin(right);
    case operator_index::cos:
        return std::cos(right);
    case operator_index::tan:
        return std::tan(right);
    case operator_index::atan:
        return std::atan(right);
    case operator_index::exp:
        return std::exp(right);
    case operator_index::abs:
        return std::abs(right);
    case operator_index::sign:
        return (right > 0) - (right < 0);
    case operator_index::sqr:
        return right * right;
    case operator_index::sqrt:
        return std::sqrt(right);
    case operator_index::log:
        return std::log(right);
    case operator_index::tgamma:
        return std::tgamma(right);
    case operator_index::lgamma:
        return std::lgamma(right);
    case operator_index::exp2:
        return std::exp2(right);
    case operator_index::expm1:
        return std::expm1(right);
    case operator_index::log10:
        return std::log10(right);
    case operator_index::log2:
        return std::log2(right);
    case operator_index::log1p:
        return std::log1p(right);
    case operator_index::cbrt:
        return std::cbrt(right);
    case operator_index::asin:
        return std::asin(right);
    case operator_index::acos:
        return std::acos(right);
    case operator_index::sinh:
        return std::sinh(right);
    case operator_index::cosh:
        return std::cosh(right);
    case operator_index::tanh:
        return std::tanh(right);
    case operator_index::asinh:
        return std::asinh(right);
    case operator_index::acosh:
        return std::acosh(right);
    case operator_index::atanh:
        return std::atanh(right);
    case operator_index::erf:
        return std::erf(right);
    case operator_index::erfc:
        return std::erfc(right);
    case operator_index::ceil:
        return std::ceil(right);
    case operator_index::floor:
        return std::floor(right);
    case operator_index::trunc:
        return std::trunc(right);
    case operator_index::round:
        return std::round(right);
    default:
        throw std::domain_error{"Error. Undefined operator <" + op + ">/."};
    }
}

template <utils::arithmetic T>
T math_expression::calc_polish_notation(const std::span<const T> input_variables) const {
    if (input_variables.size() != _variables.size()) [[unlikely]]
        throw std::domain_error{"Wrong number of variables."};
    const auto& arithmetic_operators = get_arithmetic_operators();
    std::stack<T> calculation_values;
    for (const std::string& smth : _polish_notation) {
        if (utils::is_number(smth)) {
            calculation_values.push(utils::get_number<T>(smth));
        }
        else if (_variables.count(smth) > 0) {
            calculation_values.push(input_variables[_variables.at(smth)]);
        }
        else if (arithmetic_operators.count(smth) > 0) {
            calculation_values.push(execute<T>(arithmetic_operators, calculation_values, smth));
        }
    }
    return calculation_values.top();
}

};

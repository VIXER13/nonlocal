#pragma once

#include <math_expression/math_expression.hpp>
#include <solvers/base/equation_parameters.hpp>

#include <nlohmann/json.hpp>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> read_coefficient(const nlohmann::json& config, const std::string& path) {
    static_assert(Dimension == 1 || Dimension == 2, "Supported dimensions: 1 and 2.");
    if (config.is_number())
        return config.get<T>();
    if (config.is_string()) {
        const formula::math_expression<T> parsed_formula{config.get<std::string>()};
        if constexpr (Dimension == 1) {
            if (parsed_formula.variables_count() == 1)
                return spatial_dependency<T, 1>{
                    [parsed_formula](const T argument) {
                        return parsed_formula({argument});
                    }
                };
            if (parsed_formula.variables_count() == 2)
                return solution_dependency<T, 1>{
                    [parsed_formula](const T argument, const T solution) {
                        return parsed_formula({argument, solution});
                    }
                };
        }
        if constexpr (Dimension == 2) {
            if (parsed_formula.variables_count() == 2)
                return spatial_dependency<T, 2>{
                    [parsed_formula](const std::array<T, 2>& arguments) {
                        return parsed_formula({arguments[0], arguments[1]});
                    }
                };
            if (parsed_formula.variables_count() == 3)
                return solution_dependency<T, 2>{
                    [parsed_formula](const std::array<T, 2>& arguments, const T solution) {
                        return parsed_formula({arguments[0], arguments[1], solution});
                    }
                };
        }
        throw std::domain_error{"Unsupported number of variables in coefficient \"" + path + "\"."};
    }
    throw std::domain_error{"Unsupported coefficient type in \"" + path + "\": must be a number or formula."};
}

}
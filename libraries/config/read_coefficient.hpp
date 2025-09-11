#pragma once

#include <math_expression/math_expression.hpp>
#include <solvers/base/equation_parameters.hpp>

#include <nlohmann/json.hpp>

namespace nonlocal::config {

template<std::floating_point T, size_t Dimension>
coefficient_t<T, Dimension> read_coefficient(const nlohmann::json& config, const std::string& path) {
    static_assert(Dimension == 2, "Now supported only 2D.");
    if (config.is_number())
        return config.get<T>();
    if (config.is_string()) {
        const formula::math_expression parsed_formula{config.get<std::string>()};
        if constexpr (Dimension == 2) {
            if (parsed_formula.variables_count() == 2)
                return spatial_dependency<T, Dimension>{
                    [parsed_formula](const std::array<T, 2>& arguments) { 
                        return parsed_formula({arguments[0], arguments[1]}); 
                    }
                };
            if (parsed_formula.variables_count() == 3)
                return solution_dependency<T, Dimension>{
                    [parsed_formula](const std::array<T, 2>& arguments, const T solution) { 
                        return parsed_formula({arguments[0], arguments[1], solution}); 
                    }
                };
        }
        throw std::domain_error{"Unsupported number of variables in coefficient."};
    }
    throw std::domain_error{"Unsupported coefficient type: " + path};
}

}
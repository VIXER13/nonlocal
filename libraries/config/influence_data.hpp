#pragma once

#include "config_utils.hpp"

#include <logger/logger.hpp>
#include <solvers/solver_1d/influence_functions_1d.hpp>
#include <solvers/solver_2d/influence_functions_2d.hpp>

#include <functional>
#include <variant>

namespace nonlocal::config {

enum class influence_family_t : uint8_t {
    Custom,
    Constant,
    Polynomial,
    Exponential
};

NLOHMANN_JSON_SERIALIZE_ENUM(influence_family_t, {
    {influence_family_t::Custom, nullptr},
    {influence_family_t::Constant, "constant"},
    {influence_family_t::Polynomial, "polynomial"},
    {influence_family_t::Exponential, "exponential"}
})

template<std::floating_point T>
std::variant<size_t, T> read_influence_parameter(const nlohmann::json& config, const std::string& parameter) {
    return config[parameter].is_number_integer() ? config[parameter].get<size_t>() : config[parameter].get<T>();
}

template<std::floating_point T>
std::function<T(T, T)> read_influence_1d(const nlohmann::json& config, const std::string& path, const T radius) {
    if (config.is_string())
        throw std::domain_error{"The formulaic definition of the nonlocal influence function is currently unavailable. "
                                "Please specify the function family in the " + path};
    check_required_fields(config, { "family" }, path);
    const auto family = config["family"].get<influence_family_t>();
    if (family == influence_family_t::Constant)
        return influence::constant_1d<T>{radius};
    if (family == influence_family_t::Exponential)
        return influence::normal_distribution_1d<T>{radius};
    if (family == influence_family_t::Polynomial) {
        check_required_fields(config, { "p", "q" }, path);
        const auto p = read_influence_parameter(config, "p");
        const auto q = read_influence_parameter(config, "q");
        return influence::polynomial_1d<T, 1u, 1u>(radius);
    }
    throw std::domain_error{"Unknown family of functions: " + path};
}

template<std::floating_point T>
std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> 
read_influence_1d(const nlohmann::json& config, const std::string& path, const std::array<T, 2>& radii) {
    if (config.is_string())
        throw std::domain_error{"The formulaic definition of the nonlocal influence function is currently unavailable. "
                                "Please specify the function family in the " + path};
    const std::string path_with_access = append_access_sign(path);
    check_required_fields(config, { "family" }, path);
    check_optional_fields(config, { "n" }, path);
    const auto n = config.contains("n") ? read_influence_parameter(config, "n") : 2u;
    if (config["family"].get<influence_family_t>() != influence_family_t::Constant) {
        check_required_fields(config, { "p", "q" }, path);
    }
    return {};
}

template<std::floating_point T, size_t Dimension>
struct influence_data final {
    influence_family_t family = influence_family_t::POLYNOMIAL;
    size_t n = 2;
    size_t p = 2;
    size_t q = 1;
    std::string custom;

    explicit constexpr influence_data() noexcept = default;
    explicit influence_data(const nlohmann::json& config, const std::string& path = {}) {
        const std::string path_with_access = append_access_sign(path);
        check_required_fields(config, { "family" }, path_with_access);
        family = config["family"].get<influence_family_t>();
        if (family == influence_family_t::CUSTOM) {
            custom = config["custom"].get<std::string>();
        } else {
            if constexpr (Dimension != 1) {
                check_required_fields(config, { "n" }, path_with_access);
                n = config["n"].is_number_unsigned() ?
                    config["n"].get<size_t>() :
                    config["n"].is_string() && config["n"].get<std::string>() == "infinity" ?
                    std::numeric_limits<size_t>::max() : 
                    throw std::domain_error{"Unexpected n value."};
            }
            if (family != influence_family_t::CONSTANT) {
                check_required_fields(config, { "p", "q" }, path_with_access);
                p = config["p"].get<size_t>();
                q = config["q"].get<size_t>();
            }
        }
    }

    operator nlohmann::json() const {
        if (family == influence_family_t::CUSTOM)
            return {
                {"family", "custom"},
                {"custom", custom}
            };
        nlohmann::json result = {{"family", family}};
        if constexpr (Dimension != 1) {
            if (n == std::numeric_limits<size_t>::max())
                result["n"] = "infinity";
            else
                result["n"] = n;
        }
        if (family != influence_family_t::CONSTANT) {
            result["p"] = p;
            result["q"] = q;
        }
        return result;
    }
};

// template<std::floating_point T>
// std::function<T(const std::array<T, 2>&, const std::array<T, 2>&)> get_influence(
//     const config::influence_data<T, 2>& data, const std::array<T, 2>& r) {
//     if (data.family == config::influence_family_t::CUSTOM)
//         throw std::domain_error{"Сustom influence functions are not currently supported."};
//     if (data.family == config::influence_family_t::CONSTANT) {
//         if (data.n == 1)
//             return influence::constant_2d<T, 1>{r};
//         if (data.n == 2)
//             return influence::constant_2d<T, 2>{r};
//         if (data.n == std::numeric_limits<size_t>::max())
//             return influence::constant_2d<T, std::numeric_limits<size_t>::max()>{r};
//         return influence::constant_2d_dynamic<T>{r, T(data.n)};
//     }
//     if (data.family == config::influence_family_t::POLYNOMIAL) {
//         if (data.n == 2 && data.p == 2 && data.q == 1)
//             return influence::polynomial_2d<T, 2, 1>{r};
//         if (data.n == std::numeric_limits<size_t>::max() && data.p == 1 && data.q == 1)
//             return influence::polynomial_2d<T, 1, 1, std::numeric_limits<size_t>::max()>{r};
//         return influence::polynomial_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
//     }
//     if (data.family == config::influence_family_t::EXPONENTIAL) {
//         // q == 0 temporary condition. In future refactoring will be changed 
//         if (data.n == 1 && data.p == 2 && data.q == 0)
//             return influence::exponential_2d<T, 2, 1>{r, 0.5};
//         if (data.n == 2 && data.p == 2 && data.q == 0) 
//             return influence::exponential_2d<T, 2, 2>{r, 0.5};
//         if (data.n == 2 && data.p == 3 && data.q == 0) 
//             return influence::exponential_2d<T, 3, 2>{r, 0.5};
//         if (data.n == 2 && data.p == 5 && data.q == 0) 
//             return influence::exponential_2d<T, 5, 2>{r, 0.5};
//         if (data.n == 5 && data.p == 2 && data.q == 0)
//             return influence::exponential_2d<T, 2, 5>{r, 0.5};
//         if (data.n == std::numeric_limits<size_t>::max() && data.p == 2 && data.q == 0)
//             return influence::exponential_2d<T, 2, std::numeric_limits<size_t>::max()>{r, 0.5};

//         if (data.n == 2 && data.p == 2 && data.q == 1)
//             return influence::exponential_2d<T, 2>{r, T(data.q)};
//         if (data.n == 2 && data.p == 2 && data.q == 2) // q == 2 temporary condition. In future refactoring will be changed 
//             return influence::exponential_2d<T, 2>{r, 1.5};
//         if (data.n == 2 && data.p == 2 && data.q == 3)
//             return influence::exponential_2d<T, 2>{r, T(data.q)};

//         if (data.n == 2 && data.p == 1 && data.q == 1)
//             return influence::exponential_2d<T, 1>{r, T(data.q)};
//         if (data.n == 2 && data.p == 3 && data.q == 1)
//             return influence::exponential_2d<T, 3>{r, T(data.q)};

//         return influence::exponential_2d_dynamic<T>{r, T(data.p), T(data.q), T(data.n)};
//     }
//     throw std::domain_error{"Unsupported influence family."};
// }

}
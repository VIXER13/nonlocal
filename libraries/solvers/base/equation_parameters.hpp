#pragma once

#include <constants/nonlocal_constants.hpp>
#include <metamath/types/visitor.hpp>
#include <metamath/types/vector_with_shifted_index.hpp>

#include <algorithm>
#include <array>
#include <concepts>
#include <functional>
#include <ranges>
#include <string>
#include <unordered_map>
#include <variant>

namespace nonlocal {

template<std::floating_point T, size_t Dimension>
using point = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;

template<std::floating_point T, size_t Dimension>
using spatial_dependency = std::function<T(const point<T, Dimension>&)>;

template<std::floating_point T, size_t Dimension>
using solution_dependency = std::function<T(const point<T, Dimension>&, const T)>;

template<std::floating_point T, size_t Dimension>
using coefficient_t = std::variant<T, spatial_dependency<T, Dimension>, solution_dependency<T, Dimension>>;

template<class T>
using evaluated_parameters = std::variant<T, metamath::types::vector_with_shifted_index<T>>;
constexpr size_t Constant = 0;
constexpr size_t Nonconstant = 1; // Named indices for evaluated_parameters variants

template<size_t Dimension, class T>
struct model_parameters final {
    static_assert(Dimension > 0, "Dimension must be non-zero.");
    using arg = std::conditional_t<Dimension == 1, T, std::array<T, Dimension>>;
    std::function<T(const arg&, const arg&)> influence = [](const arg&, const arg&) constexpr noexcept { return T{0}; };
    T local_weight = T{1};
};

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
struct equation_parameters final {
    model_parameters<Dimension, T> model;
    Physical<T, Args...> physical;
};

template<std::floating_point T, size_t Dimension>
bool is_constant(const coefficient_t<T, Dimension>& coefficient) {
    return std::holds_alternative<T>(coefficient);
}

template<std::floating_point T, size_t Dimension, size_t N>
bool is_constant(const std::array<coefficient_t<T, Dimension>, N>& coefficient) {
    static constexpr auto checker = [](const coefficient_t<T, Dimension>& coefficient) { return is_constant<T, Dimension>(coefficient); };
    return std::all_of(coefficient.begin(), coefficient.end(), checker);
}

template<std::floating_point T, size_t Dimension>
T evaluate(const coefficient_t<T, Dimension>& coefficient, const point<T, Dimension>& point, const T solution) {
    return std::visit(metamath::types::visitor{
        [](const T value) noexcept { return value; },
        [&point](const spatial_dependency<T, 2u>& value) { return value(point); },
        [&point, solution](const solution_dependency<T, 2u>& value) { return value(point, solution); }
    }, coefficient);
}

template<std::floating_point T, size_t Dimension, size_t N>
std::array<T, N> evaluate(const std::array<coefficient_t<T, Dimension>, N>& coefficient,
                          const point<T, Dimension>& point, const T solution) {
    std::array<T, N> result;
    for (const size_t i : std::ranges::iota_view{0zu, N})
        result[i] = evaluate<T, Dimension>(coefficient[i], point, solution);
    return result;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::vector<model_parameters<Dimension, T>> get_models(const std::vector<equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::vector<model_parameters<Dimension, T>> models(parameters.size());
    std::transform(parameters.cbegin(), parameters.cend(), models.begin(),
                  [](const equation_parameters<Dimension, T, Physical, Args...>& parameter) { return parameter.model; });
    return models;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::unordered_map<std::string, model_parameters<Dimension, T>> get_models(
    const std::unordered_map<std::string, equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::unordered_map<std::string, model_parameters<Dimension, T>> models;
    for(const auto& [name, parameter] : parameters)
        models[name] = parameter.model;
    return models;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::vector<Physical<T, Args...>> get_physical_parameters(const std::vector<equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::vector<Physical<T, Args...>> physical_parameters(parameters.size());
    std::transform(parameters.cbegin(), parameters.cend(), physical_parameters.begin(),
                  [](const equation_parameters<Dimension, T, Physical, Args...>& parameter) { return parameter.physical; });
    return physical_parameters;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::unordered_map<std::string, Physical<T, Args...>> get_physical_parameters(
    const std::unordered_map<std::string, equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::unordered_map<std::string, Physical<T, Args...>> physical_parameters;
    for(const auto& [name, parameter] : parameters)
        physical_parameters[name] = parameter.physical;
    return physical_parameters;
}

template<size_t Dimension, class T>
std::vector<theory_t> theories_types(const std::vector<model_parameters<Dimension, T>>& models) {
    std::vector<theory_t> theories(models.size());
    std::transform(models.cbegin(), models.cend(), theories.begin(), 
                   [](const model_parameters<Dimension, T>& model) { return theory_type(model.local_weight); } );
    return theories;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::vector<theory_t> theories_types(const std::vector<equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::vector<theory_t> theories(parameters.size());
    std::transform(parameters.cbegin(), parameters.cend(), theories.begin(), 
                   [](const equation_parameters<Dimension, T, Physical, Args...>& parameter)
                   { return theory_type(parameter.model.local_weight); } );
    return theories;
}

template<size_t Dimension, class T, template<class, auto...> class Physical, auto... Args>
std::unordered_map<std::string, theory_t> theories_types(const std::unordered_map<std::string, equation_parameters<Dimension, T, Physical, Args...>>& parameters) {
    std::unordered_map<std::string, theory_t> theories;
    for(const auto& [name, parameter] : parameters)
        theories[name] = theory_type(parameter.model.local_weight);
    return theories;
}

}
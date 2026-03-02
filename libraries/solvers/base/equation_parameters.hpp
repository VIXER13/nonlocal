#pragma once

#include "coefficient.hpp"

#include <constants/nonlocal_constants.hpp>
#include <metamath/types/vector_with_shifted_index.hpp>

#include <string>
#include <unordered_map>

namespace nonlocal {

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
#pragma once

#include <metamath/metamath.hpp>
#include "../../equation_parameters.hpp"

#include <string>
#include <unordered_map>
#include <variant>

namespace nonlocal::thermal {

template<std::floating_point T>
using spatial_dependency = std::function<T(const std::array<T, 2>&)>;

template<std::floating_point T>
using solution_dependency = std::function<T(const T, const std::array<T, 2>&)>;

template<std::floating_point T>
using coefficient_t = std::variant<T, spatial_dependency<T>, solution_dependency<T>>;

template<std::floating_point T>
using isotropic_conductivity_t = coefficient_t<T>;

template<std::floating_point T>
using orthotropic_conductivity_t = std::array<coefficient_t<T>, 2>;

template<std::floating_point T>
using anisotropic_conductivity_t = metamath::types::square_matrix<coefficient_t<T>, 2>;

template<std::floating_point T>
using conductivity_t = std::variant<
    isotropic_conductivity_t<T>,
    orthotropic_conductivity_t<T>,
    anisotropic_conductivity_t<T>
>;

template<std::floating_point T>
struct parameter_2d final {
    conductivity_t<T> conductivity = T{1};
    T capacity = T{1};
    T density = T{1};
    T relaxation_time = T{0};
};

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

template<class... Ts>
struct visitor : Ts... {
    using Ts::operator()...;
};

template<class... Ts>
visitor(Ts...) -> visitor<Ts...>;

}
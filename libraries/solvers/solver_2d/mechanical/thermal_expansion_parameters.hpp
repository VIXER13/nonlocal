#pragma once

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::mechanical {

// single constant
template<class T>
using isotropic_thermal_expansion_t = T;
// a[0]  0
//  0   a[1]
template<class T>
using orthotropic_thermal_expansion_t = std::array<T, 2>;
// a[0] a[2]
// a[2] a[1]
template<class T>
using anisotropic_thermal_expansion_t = std::array<T, 3>;
template<class T>
using thermal_expansion_t = std::variant<
    isotropic_thermal_expansion_t<T>,
    orthotropic_thermal_expansion_t<T>,
    anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
using raw_isotropic_thermal_expansion_t = isotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_orthotropic_thermal_expansion_t = orthotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_anisotropic_thermal_expansion_t = anisotropic_thermal_expansion_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_thermal_expansion_t = std::variant<
    std::monostate,
    raw_isotropic_thermal_expansion_t<T>,
    raw_orthotropic_thermal_expansion_t<T>,
    raw_anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
struct isotropic_thermal_strain final {
    metamath::types::vector_with_shifted_index<T> strain;

    std::array<T, 3> operator[](const size_t qshift) const {
        const T result = strain[qshift];
        return {result, result, T{0}};
    }
};

template<std::floating_point T>
struct orthotropic_thermal_strain final {
    metamath::types::vector_with_shifted_index<std::array<T, 2>> strain;

    std::array<T, 3> operator[](const size_t qshift) const {
        const auto& result = strain[qshift];
        return {result[0], result[1], T{0}};
    }
};

template<std::floating_point T>
struct anisotropic_thermal_strain final {
    metamath::types::vector_with_shifted_index<std::array<T, 3>> strain;

    std::array<T, 3> operator[](const size_t qshift) const {
        return strain[qshift];
    }
};

template<std::floating_point T>
struct orthotropic_constant_thermal_strain final {
    metamath::types::vector_with_shifted_index<T> delta_temperature;
    std::array<T, 2> thermal_expansion;

    std::array<T, 3> operator[](const size_t qshift) const {
        using namespace metamath::functions;
        const auto result = thermal_expansion * delta_temperature[qshift];
        return {result[0], result[1], T{0}};
    }
};

template<std::floating_point T>
struct anisotropic_constant_thermal_strain final {
    metamath::types::vector_with_shifted_index<T> delta_temperature;
    std::array<T, 3> thermal_expansion;

    std::array<T, 3> operator[](const size_t qshift) const {
        using namespace metamath::functions;
        return thermal_expansion * delta_temperature[qshift];
    }
};

template<std::floating_point T>
using evaluated_thermal_strain_t = std::variant<
    std::monostate,
    isotropic_thermal_strain<T>,
    orthotropic_thermal_strain<T>,
    anisotropic_thermal_strain<T>,
    orthotropic_constant_thermal_strain<T>,
    anisotropic_constant_thermal_strain<T>
>;

}
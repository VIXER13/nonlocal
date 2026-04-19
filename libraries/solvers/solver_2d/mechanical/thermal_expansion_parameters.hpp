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
    raw_isotropic_thermal_expansion_t<T>,
    raw_orthotropic_thermal_expansion_t<T>,
    raw_anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
using evaluated_isotropic_thermal_expansion_t = evaluated_parameters<isotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_orthotropic_thermal_expansion_t = evaluated_parameters<orthotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_anisotropic_thermal_expansion_t = evaluated_parameters<anisotropic_thermal_expansion_t<T>>;
template<std::floating_point T>
using evaluated_thermal_expansion_t = std::variant<
    evaluated_isotropic_thermal_expansion_t<T>,
    evaluated_orthotropic_thermal_expansion_t<T>,
    evaluated_anisotropic_thermal_expansion_t<T>
>;

template<std::floating_point T>
using isotropic_thermal_strain_t = metamath::types::vector_with_shifted_index<T>;
template<std::floating_point T>
using orthotropic_thermal_strain_t = metamath::types::vector_with_shifted_index<std::array<T, 2>>;
template<std::floating_point T>
using anisotropic_thermal_strain_t = metamath::types::vector_with_shifted_index<std::array<T, 3>>;
template<std::floating_point T>
using orthotropic_constant_thermal_strain_t = std::pair<metamath::types::vector_with_shifted_index<T>, std::array<T, 2>>;
template<std::floating_point T>
using anisotropic_constant_thermal_strain_t = std::pair<metamath::types::vector_with_shifted_index<T>, std::array<T, 3>>;
template<std::floating_point T>
using evaluated_thermal_strain_t = std::variant<
    isotropic_thermal_strain_t<T>,
    orthotropic_thermal_strain_t<T>,
    anisotropic_thermal_strain_t<T>,
    orthotropic_constant_thermal_strain_t<T>,
    anisotropic_constant_thermal_strain_t<T>
>;

template<std::floating_point T>
struct evaluated_thermal_strain final {
    evaluated_thermal_strain_t<T> strain = {};

    std::array<T, 3> operator[](const size_t qshift) const {
        return std::visit(metamath::types::visitor{
            [qshift](const isotropic_thermal_strain_t<T>& strain) {
                const T result = strain[qshift];
                return std::array{result, result, T{0}};
            },
            [qshift](const orthotropic_thermal_strain_t<T>& strain) {
                const auto& result = strain[qshift];
                return std::array{result[XX], result[YY], T{0}};
            },
            [qshift](const anisotropic_thermal_strain_t<T>& strain) {
                return strain[qshift];
            },
            [qshift](const orthotropic_constant_thermal_strain_t<T>& strain) {
                using namespace metamath::functions;
                const auto& [delta_temperature, thermal_expansion] = strain;
                const auto result = thermal_expansion * delta_temperature[qshift];
                return std::array{result[XX], result[YY], T{0}};
            },
            [qshift](const anisotropic_constant_thermal_strain_t<T>& strain) {
                using namespace metamath::functions;
                const auto& [delta_temperature, thermal_expansion] = strain;
                return thermal_expansion * delta_temperature[qshift];
            }
        }, strain);
    }
};

}
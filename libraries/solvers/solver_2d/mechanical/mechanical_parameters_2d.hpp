#pragma once

#include <metamath/functions/power.hpp>
#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::mechanical {

// d[0] d[1]  0
// d[1] d[0]  0
//  0    0   d[2]
template<class T>
using isotropic_hook_matrix_t = std::array<T, 3>;
// d[0] d[1]  0
// d[1] d[2]  0
//  0    0   d[3]
template<class T>
using orthotropic_hook_matrix_t = std::array<T, 4>;
// d[0] d[1] d[2]
// d[1] d[3] d[4]
// d[2] d[4] d[5]
template<class T>
using anisotropic_hook_matrix_t = std::array<T, 6>;
template<class T>
using hooke_matrix_t = std::variant<
    isotropic_hook_matrix_t<T>,
    orthotropic_hook_matrix_t<T>,
    anisotropic_hook_matrix_t<T>
>;

template<std::floating_point T>
using raw_isotropic_hook_matrix_t = isotropic_hook_matrix_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_orthotropic_hook_matrix_t = orthotropic_hook_matrix_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_anisotropic_hook_matrix_t = anisotropic_hook_matrix_t<coefficient_t<T, 2>>;
template<std::floating_point T>
using raw_hook_matrix_t = std::variant<
    raw_isotropic_hook_matrix_t<T>,
    raw_orthotropic_hook_matrix_t<T>,
    raw_anisotropic_hook_matrix_t<T>
>;

template<std::floating_point T>
using evaluated_isotropic_hook_matrix_t = evaluated_parameters<isotropic_hook_matrix_t<T>>;
template<std::floating_point T>
using evaluated_orthotropic_hook_matrix_t = evaluated_parameters<orthotropic_hook_matrix_t<T>>;
template<std::floating_point T>
using evaluated_anisotropic_hook_matrix_t = evaluated_parameters<anisotropic_hook_matrix_t<T>>;
template<std::floating_point T>
using evaluated_hook_matrix_t = std::variant<
    evaluated_isotropic_hook_matrix_t<T>,
    evaluated_orthotropic_hook_matrix_t<T>,
    evaluated_anisotropic_hook_matrix_t<T>
>;

template<std::floating_point T>
using evaluated_hook_matrices_2d = std::unordered_map<std::string, equation_parameters<2, T, evaluated_hook_matrix_t>>;

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
using evaluated_thermal_expansion_2d = std::unordered_map<std::string, equation_parameters<2, T, evaluated_thermal_expansion_t>>;

template<std::floating_point T>
struct isotropic_elastic_parameters final {
    coefficient_t<T, 2> young_modulus = T{210.};
    coefficient_t<T, 2> poissons_ratio = T{0.3};
    raw_isotropic_thermal_expansion_t<T> thermal_expansion = T{13e-6};

    isotropic_hook_matrix_t<T> hooke(const std::array<T, 2>& x) const {
        const T E = evaluate<T, 2zu>(young_modulus, x, {});
        const T nu = evaluate<T, 2zu>(poissons_ratio, x, {});
        const T value = E / (T{1} - nu * nu);
        return { value, nu * value, 0.5 * value * (1 - nu) };
    }
};

template<std::floating_point T>
struct orthotropic_elastic_parameters final {
    std::array<coefficient_t<T, 2>, 2> young_modulus = {T{210.}, T{210.}};
    std::array<coefficient_t<T, 2>, 2> poissons_ratio = {T{0.3}, T{0.3}};
    coefficient_t<T, 2> shear_modulus = T{210.} / T{13};
    raw_orthotropic_thermal_expansion_t<T> thermal_expansion = {T{13e-6}, T{13e-6}};

    orthotropic_hook_matrix_t<T> hooke(const std::array<T, 2>& x) const {
        const auto E = evaluate<T, 2zu>(young_modulus, x, {});
        const auto nu = evaluate<T, 2zu>(poissons_ratio, x, {});
        const T G = evaluate<T, 2zu>(shear_modulus, x, {});
        const T value = T{1} / (T{1} - nu[0] * nu[1]);
        return {E[X] * value, E[X] * nu[1] * value, E[Y] * value, G};
    }
};

template<std::floating_point T>
struct anisotropic_elastic_parameters final {
    orthotropic_elastic_parameters<T> main_parameters;
    coefficient_t<T, 2> angle = T{0};
    raw_anisotropic_thermal_expansion_t<T> thermal_expansion = {T{13e-6}, T{13e-6}, T{0}};

    static anisotropic_hook_matrix_t<T> rotate(const orthotropic_hook_matrix_t<T>& matrix, const T angle) noexcept {
        const T sin = std::sin(angle);
        const T cos = std::cos(angle);
        return {0, 1, 2, 3, 4, 5};
    }

    anisotropic_hook_matrix_t<T> hooke(const std::array<T, 2>& x) const {
        return rotate(main_parameters.hooke(x), evaluate<T, 2zu>(angle, x, {}));
    }
    
    operator raw_anisotropic_thermal_expansion_t<T>&() const {
        return main_parameters;
    }
};

template<std::floating_point T>
using elastic_parameters_t = std::variant<
    isotropic_elastic_parameters<T>,
    orthotropic_elastic_parameters<T>,
    anisotropic_elastic_parameters<T>
>;

template<std::floating_point T>
using elastic_parameters = std::unordered_map<std::string, equation_parameters<2, T, elastic_parameters_t>>;

}
#pragma once

#include <solvers/base/equation_parameters.hpp>

namespace nonlocal::solver_2d::mechanical {

namespace isotropic_indices   { enum : uint8_t {_11, _12, _66}; }
namespace orthotropic_indices { enum : uint8_t {_11, _12, _22, _66}; }
namespace anisotropic_indices { enum : uint8_t {_11, _12, _16, _22, _26, _66}; }

// d[_11] d[_12]   0
// d[_12] d[_11]   0
//   0      0    d[_66]
template<class T>
using isotropic_hook_matrix_t = std::array<T, 3>;
// d[_11] d[_12]   0
// d[_12] d[_22]   0
//   0      0    d[_66]
template<class T>
using orthotropic_hook_matrix_t = std::array<T, 4>;
// d[_11] d[_12] d[_16]
// d[_12] d[_22] d[_26]
// d[_16] d[_26] d[_66]
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
struct isotropic_elastic_parameters final {
    coefficient_t<T, 2> young_modulus = T{210.};
    coefficient_t<T, 2> poissons_ratio = T{0.3};

    bool is_constant() const {
        return nonlocal::is_constant(young_modulus) && nonlocal::is_constant(poissons_ratio);
    }

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

    bool is_constant() const {
        return nonlocal::is_constant(young_modulus) && 
               nonlocal::is_constant(poissons_ratio) &&
               nonlocal::is_constant(shear_modulus);
    }

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

    bool is_constant() const {
        return main_parameters.is_constant() && nonlocal::is_constant(angle);
    }

    static anisotropic_hook_matrix_t<T> rotate(const orthotropic_hook_matrix_t<T>& matrix, const T angle) noexcept {
        const T sin = std::sin(angle);
        const T cos = std::cos(angle);
        const T sin2 = sin * sin;
        const T cos2 = cos * cos;
        const T sin4 = sin2 * sin2;
        const T cos4 = cos2 * cos2;
        const T sin2cos2 = sin2 * cos2;
        const T sin3cos = sin2 * sin * cos;
        const T sincos3 = sin * cos * cos2;
        using namespace orthotropic_indices;
        return {
            matrix[_11] * cos4 + matrix[_22] * sin4 + (2 * matrix[_12] + 4 * matrix[_66]) * sin2cos2,
            (matrix[_11] + matrix[_22] - 4 * matrix[_66]) * sin2cos2 + matrix[_12] * (cos4 + sin4),
            (matrix[_11] - matrix[_12] - 2 * matrix[_66]) * sincos3 - (matrix[_22] - matrix[_12] - 2 * matrix[_66]) * sin3cos,
            matrix[_11] * sin4 + matrix[_22] * cos4 + (2 * matrix[_12] + 4 * matrix[_66]) * sin2cos2,
            (matrix[_11] - matrix[_12] - 2 * matrix[_66]) * sin3cos - (matrix[_22] - matrix[_12] - 2 * matrix[_66]) * sincos3,
            (matrix[_11] + matrix[_22] - 2 * (matrix[_12] + matrix[_66])) * sin2cos2 + matrix[_66] * (cos4 + sin4)
        };
    }

    anisotropic_hook_matrix_t<T> hooke(const std::array<T, 2>& x) const {
        return rotate(main_parameters.hooke(x), evaluate<T, 2zu>(angle, x, {}));
    }
};

template<std::floating_point T>
using elastic_parameters_t = std::variant<
    isotropic_elastic_parameters<T>,
    orthotropic_elastic_parameters<T>,
    anisotropic_elastic_parameters<T>
>;

template<std::floating_point T>
bool is_constant(const elastic_parameters_t<T>& elastic) {
    return std::visit([](const auto& elastic) { return elastic.is_constant(); }, elastic);
}

}
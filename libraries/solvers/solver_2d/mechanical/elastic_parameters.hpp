#pragma once

#include <solvers/base/equation_parameters.hpp>

#include <cmath>

namespace nonlocal::solver_2d::mechanical {

namespace isotropic_indices   { enum : uint8_t {_11, _12, _66}; }
namespace orthotropic_indices { enum : uint8_t {_11, _12, _66, _22}; }
namespace anisotropic_indices { enum : uint8_t {_11, _12, _66, _22, _16, _26}; }

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
        return { value, nu * value, T{0.5} * value * (T{1} - nu) };
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
        return {E[X] * value, E[X] * nu[1] * value, G, E[Y] * value};
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
        using namespace orthotropic_indices;
        const T angle_2 = angle + angle;
        const T sin_2 = std::sin(angle_2);
        const T cos_2 = std::cos(angle_2);
        const T cos_4 = T{0.125} * (cos_2 * cos_2 - sin_2 * sin_2);
        const T stiff_diff = matrix[_11] - matrix[_22];
        const T stiff_core = matrix[_11] + matrix[_22] - 2 * matrix[_12] - 4 * matrix[_66];
        const T main_part_plus  = stiff_core * (cos_4 + T{0.125});
        const T main_part_minus = stiff_core * (T{0.125} - cos_4);
        const T cos_diff = stiff_diff * (cos_2 + T{0.5});
        const T minor_part = stiff_core * cos_2;
        const T sin_2_025 = T{0.25} * sin_2;
        return {
            main_part_plus + T{0.5} * (matrix[_12] + matrix[_22] + cos_diff) + matrix[_66],
            main_part_minus + matrix[_12],
            main_part_minus + matrix[_66],
            main_part_plus + T{0.5} * (matrix[_12] + matrix[_11] - cos_diff) + matrix[_66],
            sin_2_025 * (stiff_diff + minor_part),
            sin_2_025 * (stiff_diff - minor_part)
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
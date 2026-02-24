#pragma once

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
using hook_matrix_t = std::variant<
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
using evaluated_hook_matrix_2d = std::unordered_map<std::string, equation_parameters<2, T, evaluated_hook_matrix_t>>;

enum class plane_t : bool {
    STRESS,
    STRAIN
}; 

// arr[0] arr[1]   0
// arr[1] arr[0]   0
//   0      0    arr[2]
template<class T>
using hooke_matrix = std::array<T, 4>;

template<class T>
struct parameter_2d final {
    std::array<T, 2> youngs_modulus = {210, 210};
    std::array<T, 2> poissons_ratio = {0.3, 0.3};
    T shear_modulus = 80;
    T thermal_expansion = 13e-6;

    material_t material = material_t::ISOTROPIC;

    constexpr T E(const plane_t plane, const std::size_t i) const noexcept;
    constexpr T nu(const plane_t plane, const std::size_t i) const noexcept;
    constexpr T G() const noexcept;
    constexpr hooke_matrix<T> hooke(const plane_t plane) const noexcept;
};

template<class T>
using parameters_2d = std::unordered_map<std::string, equation_parameters<2, T, parameter_2d>>;

template<class T>
struct mechanical_parameters_2d final {
    parameters_2d<T> materials;
    std::vector<T> delta_temperature; // if empty, then thermal expansions are not taken into account
    plane_t plane = plane_t::STRESS;
};

template<class T>
constexpr T parameter_2d<T>::E(const plane_t plane, const std::size_t i) const noexcept {
    return plane == plane_t::STRESS ? youngs_modulus[i] : youngs_modulus[i] / (T{1} - poissons_ratio[i] * poissons_ratio[i]);
}

template<class T>
constexpr T parameter_2d<T>::nu(const plane_t plane, const std::size_t i) const noexcept {
    return plane == plane_t::STRESS ? poissons_ratio[i] : poissons_ratio[i] / (T{1} - poissons_ratio[i]);
}

template<class T>
constexpr T parameter_2d<T>::G() const noexcept {
    return shear_modulus;
}

template<class T>
constexpr hooke_matrix<T> parameter_2d<T>::hooke(const plane_t plane) const noexcept {
    const T Ex = this->E(plane, 0);
    const T Ey = this->E(plane, 1);
    const T nu_xy = this->nu(plane, 0);
    const T nu_yx = this->nu(plane, 1);
    const T Gxy = material == material_t::ISOTROPIC ? T{0.5} * Ex / (1 + nu_xy) : this->G();
    const T div = T{1} / (T{1} - nu_xy*nu_yx);
    // Ey * nu_xy == Ex * nu_yx
    return { Ex * div, Ex * nu_yx * div,
             Ey * div, Gxy};
}

}
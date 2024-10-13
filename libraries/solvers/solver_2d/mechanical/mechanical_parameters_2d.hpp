#ifndef NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP
#define NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP

#include "../../equation_parameters.hpp"

#include <array>
#include <string>
#include <unordered_map>

namespace nonlocal::mechanical {

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
    matrix_init_t matrix_init = matrix_init_t::X_dominant;

    constexpr T E(const plane_t plane, const std::size_t i) const noexcept;
    constexpr T nu(const plane_t plane, const std::size_t i) const noexcept;
    constexpr T G() const noexcept;
    constexpr hooke_matrix<T> hooke(const plane_t plane) const noexcept;
    constexpr T side_diagonal(const T Ex, const T Ey, const T nuyx, const T nuxy, const T div) const noexcept;
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
    return material == material_t::ISOTROPIC ? T{0.5} * youngs_modulus[0] / (1 + poissons_ratio[0]) : shear_modulus;
}

template<class T>
constexpr T parameter_2d<T>::side_diagonal(const T Ex, const T Ey, const T nuyx, const T nuxy, const T div) const noexcept {
    return matrix_init == matrix_init_t::X_dominant ? Ex * nuyx * div  : Ey * nuxy * div;
}

template<class T>
constexpr hooke_matrix<T> parameter_2d<T>::hooke(const plane_t plane) const noexcept {
    const T Ex = this->E(plane, 0);
    const T Ey = this->E(plane, 1);
    const T nuxy = this->nu(plane, 0);
    const T nuyx = this->nu(plane, 1);
    const T Gxy = this->G();
    const T div = T{1} / (T{1} - nuxy*nuyx);
    // Ey * nuxy = Ex * nuyx
    return { Ex * div, side_diagonal(Ex, Ey, nuyx, nuxy, div),
             Ey * div, Gxy};
}

}

#endif
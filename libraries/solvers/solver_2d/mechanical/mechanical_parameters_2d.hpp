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
using hooke_matrix = std::array<T, 3>;

template<class T>
struct parameter_2d final {
    T young_modulus = 210;
    T poissons_ratio = 0.3;
    T thermal_expansion = 13e-6;

    constexpr T E(const plane_t plane) const noexcept;
    constexpr T nu(const plane_t plane) const noexcept;
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
constexpr T parameter_2d<T>::E(const plane_t plane) const noexcept {
    return plane == plane_t::STRESS ? young_modulus : young_modulus / (T{1} - poissons_ratio * poissons_ratio);
}

template<class T>
constexpr T parameter_2d<T>::nu(const plane_t plane) const noexcept {
    return plane == plane_t::STRESS ? poissons_ratio : poissons_ratio / (T{1} - poissons_ratio);
}

template<class T>
constexpr hooke_matrix<T> parameter_2d<T>::hooke(const plane_t plane) const noexcept {
    const T E = this->E(plane);
    const T nu = this->nu(plane);
    return {     E / (T{1} - nu*nu),
            nu * E / (T{1} - nu*nu),
        T{0.5} * E / (T{1} + nu) };
}

}

#endif
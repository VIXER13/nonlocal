#ifndef NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP
#define NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP

#include "functions/power.hpp"

#include <array>
#include <vector>

namespace nonlocal::mechanical {

enum class task_2d_t : bool { PLANE_STRESS, PLANE_STRAIN };

template<class T>
struct equation_parameters final {
    task_2d_t task = task_2d_t::PLANE_STRESS;
    T poissons_ratio = 0.3;
    T young_modulus = 210;
    T thermal_expansion = 13e-6;
    bool is_thermoelasticity = false;
    std::vector<T> delta_temperature;

    constexpr T E() const noexcept {
        if (task == task_2d_t::PLANE_STRAIN)
            return young_modulus / (T{1} - metamath::functions::power<2>(poissons_ratio));
        return young_modulus;
    }

    constexpr T nu() const noexcept {
        if (task == task_2d_t::PLANE_STRAIN)
            return poissons_ratio / (T{1} - poissons_ratio);
        return poissons_ratio;
    }
};

// Матрица Гука, которая имеет следующий портрет:
// arr[0] arr[1]   0
// arr[1] arr[0]   0
//   0      0    arr[2]
template<class T>
constexpr std::array<T, 3> hooke_matrix(const equation_parameters<T>& parameters) noexcept {
    const T nu = parameters.nu(),
            E  = parameters.E();
    return {       E / (1 - nu*nu),
              nu * E / (1 - nu*nu),
             0.5 * E / (1 + nu) };
}

}

#endif
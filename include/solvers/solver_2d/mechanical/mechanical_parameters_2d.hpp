#ifndef NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP
#define NONLOCAL_MECHANICAL_PARAMETERS_2D_HPP

#include <array>
#include <vector>

namespace nonlocal::mechanical {

enum class calc_t : bool { PLANE_STRESS, PLANE_STRAIN };

template<class T>
struct equation_parameters final {
    calc_t type = calc_t::PLANE_STRESS;
    T nu    = 0, // Коэффициент Пуассона
    E     = 0, // Модуль Юнга
    alpha = 0, // Коэффициент линейного расширения
    p1    = 0, // Весовой параметр модели
    r     = 0; // Длины полуосей области нелокального влияния
    bool thermoelasticity = false; // Учитывать температурные деформации
    std::vector<T> delta_temperature; // Разница температур: T - T0

    T poisson() const noexcept { return type == calc_t::PLANE_STRESS ? nu : nu / (1 - nu); }
    T young  () const noexcept { return type == calc_t::PLANE_STRESS ? E  : E  / (1 - nu * nu); }
};

// Матрица Гука, которая имеет следующий портрет:
// arr[0] arr[1]   0
// arr[1] arr[0]   0
//   0      0    arr[2]
template<class T>
constexpr std::array<T, 3> hooke_matrix(const equation_parameters<T>& parameters) noexcept {
    const T nu = parameters.poisson(),
            E  = parameters.young();
    return {       E / (1 - nu*nu),
                   nu * E / (1 - nu*nu),
                   0.5 * E / (1 + nu) };
}

}

#endif
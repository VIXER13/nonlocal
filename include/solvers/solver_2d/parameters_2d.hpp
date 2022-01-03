#ifndef NONLOCAL_PARAMETERS_2D_HPP
#define NONLOCAL_PARAMETERS_2D_HPP

#include "../solvers_constants.hpp"
#include <array>
#include <unordered_map>
#include <filesystem>

namespace nonlocal {

namespace thermal {

template<class T, material_t Material>
struct equation_parameters final {
    static_assert(Material == material_t::ISOTROPIC ||
                  Material == material_t::ORTHOTROPIC,
                  "Only isotropic and orthotropic materials are supported");

    T c        = T{1}; // Коэффициент теплоёмкости
    T rho      = T{1}; // Плотность материала
    T integral = T{0}; // Интеграл от решения (для задачи Неймана)
    std::conditional_t<Material == material_t::ISOTROPIC, T, std::array<T, 2>> lambda; // Коэффициент теплопроводности
    std::unordered_map<std::string, T> alpha; // Коэффициент теплоотдачи

    equation_parameters() noexcept {
        if constexpr(Material == material_t::ISOTROPIC)
            lambda = T{1};
        if constexpr(Material == material_t::ORTHOTROPIC)
            lambda.fill(T{1});
    }
};

}

}

#endif
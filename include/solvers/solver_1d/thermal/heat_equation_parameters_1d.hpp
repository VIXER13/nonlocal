#ifndef HEAT_EQUATION_PARAMETERS_1D_HPP
#define HEAT_EQUATION_PARAMETERS_1D_HPP

#include <array>

namespace nonlocal::thermal {

template<class T>
struct heat_equation_parameters_1d final {
    T lambda   = T{1};                     // Коэффициент теплопроводности
    T c        = T{1};                     // Коэффициент теплоёмкости
    T rho      = T{1};                     // Плотность материала
    T integral = T{0};                     // Интеграл от решения (для задачи Неймана)
    std::array<T, 2> alpha = {T{0}, T{0}}; // Коэффициент теплоотдачи на левой и правой границах
};

}

#endif
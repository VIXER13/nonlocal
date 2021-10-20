#ifndef PARAMETERS_1D_HPP
#define PARAMETERS_1D_HPP

#include <array>
#include <filesystem>

namespace nonlocal {

template<class T>
struct nonlocal_parameters_1d final {
    T p1 = T{1}; // Весовой параметр модели
    T r  = T{0}; // Радиус нелокального влияния
};

template<class T>
struct nonstationary_solver_parameters_1d final {
    std::filesystem::path save_path = std::filesystem::current_path();
    std::array<T, 2> time_interval = {0, 1};
    uintmax_t steps = 100;
    uintmax_t save_freq = 1; // Частота сохранения
    bool save_csv    = true; // Сохранять .csv файлы в формате (x, T)
    bool calc_energy = true; // Вычислять энергия при сохранении, иногда полезно для контроля расчёта
};

namespace thermal {

template<class T>
struct heat_equation_parameters_1d final {
    T lambda   = T{1};                     // Коэффициент теплопроводности
    T c        = T{1};                     // Коэффициент теплоёмкости
    T rho      = T{1};                     // Плотность материала
    T integral = T{0};                     // Интеграл от решения (для задачи Неймана)
    std::array<T, 2> alpha = {T{0}, T{0}}; // Коэффициент теплоотдачи на левой и правой границах
};

}

}

#endif
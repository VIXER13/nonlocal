#ifndef HEAT_EQUATION_PARAMETERS_HPP
#define HEAT_EQUATION_PARAMETERS_HPP

#include <array>

namespace nonlocal::heat {

enum class material_t : bool { ISOTROPIC, ORTHOTROPIC };

template<class T, material_t Material>
struct equation_parameters final {
    std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1> lambda, // Коэффициент теплопроводности
                                                               r;      // Длины полуосей области нелокального влияния
    T rho      = T{1}, // Плотность
      c        = T{1}, // Теплоёмкость
      p1       = T{1}, // Весовой параметр модели
      integral = T{0}; // Значение интеграла для задачи Неймана

    equation_parameters() noexcept {
        lambda.fill(T{1});
        r.fill(T{0});
    }
};

template<class T>
struct solver_parameters final {
    std::string save_path; // Путь куда сохранять данные
    std::array<T, 2> time_interval = {0, 1};
    uintmax_t steps = 100,
              save_freq = 1; // Частота сохранения
    bool save_vtk = true,    // Сохранять .vtk файлы
         save_csv = true,    // Сохранять .csv файлы в формате (x1, x2, T)
         calc_energy = true; // Вычислять энергия при сохранении, иногда полезно для контроля расчёта
};

}

#endif
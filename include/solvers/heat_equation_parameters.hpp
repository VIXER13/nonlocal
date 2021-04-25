#ifndef HEAT_EQUATION_PARAMETERS_HPP
#define HEAT_EQUATION_PARAMETERS_HPP

#include <array>

namespace nonlocal::heat {

enum class calc_type : bool { ISOTROPIC, ORTHOTROPIC };

template<class T, calc_type Calc_Type>
struct equation_parameters final {
    std::array<T, Calc_Type == calc_type::ORTHOTROPIC ? 2 : 1> lambda;
    T rho = 1,
      c   = 1;

    equation_parameters() {
        lambda.fill(T{1});
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
#ifndef FINITE_ELEMENT_2D_BASIS_QUADRATIC_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_BASIS_QUADRATIC_SERENDIPITY_HPP

// Базисные функции элемента взяты из статьи:
// Астионенко И.А., Гучек П.И., Литвиненко Е.И., Хомченко А.Н. Применениче альтренативных серендиповых моделей при решении задачи о кручении призматических стержней.
// Статьи с построением данного базиса я не нашёл.

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class quadratic_serendipity : public geometry_2d<T, rectangle_element_geometry> {
public:
    T get_parameter() const noexcept { return _p; }
    void set_parameter(const T p) noexcept { _p = p; }

protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;
    
    explicit quadratic_serendipity() = default;
    ~quadratic_serendipity() override = default;

    // В серендиповой аппроксимации высших порядков возникает проблема с негативизмом стандартного базиса в угловых узлах.
    // Для этого вводится специальный параметр p, который позволяет её избежать.
    // В сущности p является значением интеграла по области элемента от угловой функции. Значение интегралов от промежуточных функций есть 1-p.
    //T _p = T{-1} / T{3}; // Значение по умолчанию даёт нам классический вариант квадратичных серендиповых элементов.
    T _p = T{2} / T{9}; // Данное значение обеспечивает минимальный след матрицы, что ускоряет сходимость решателей СЛАУ
    static constexpr symdiff::variable<2> p{}; // В силу особенностей вычисления выражений, дополнительному параметру необходимо завести дополнительную переменную.

    // Нумерация узлов на квадратичном серендиповом элементе: 6---5---4
    //                                                        |       |
    //                                                        7       3
    //                                                        |       |
    //                                                        0---1---2
    static constexpr std::array<std::array<T, 2>, 8> nodes = { T{-1}, T{-1},
                                                               T{ 0}, T{-1},
                                                               T{ 1}, T{-1},
                                                               T{ 1}, T{ 0},
                                                               T{ 1}, T{ 1},
                                                               T{ 0}, T{ 1},
                                                               T{-1}, T{ 1},
                                                               T{-1}, T{ 0} };

    // Базисные функции в локальной системе координат имеют вид:
    // N_i = 0.0625 (1 + xi_i xi)(1 + eta_i eta)[(36p-1)(1 - xi_i xi - eta_i eta) + (36p+3)xi_i xi eta_i eta], xi_i = +-1, eta_i = +-1, i = 0,2,4,6,
    // N_i = 0.0625 (1 -  xi^2)(1 + eta_i eta)[(5-36p) + (36p+3)eta_i eta], eta_i = +-1, i = 1,5,
    // N_i = 0.0625 (1 - eta^2)(1 +  xi_i  xi)[(5-36p) + (36p+3) xi_i  xi],  xi_i = +-1, i = 3,7.
    static constexpr auto basis = std::make_tuple(
         (T{1}-xi)      * (T{1}-eta) * ((T{0.5625}*p-T{0.0625})*(T{1}+xi+eta) + (T{0.5625}*p+T{0.1875})*xi*eta),
        -(T{1}-xi*xi)   * (T{1}-eta) * ((T{0.5625}*p-T{0.3125})               + (T{0.5625}*p+T{0.1875})*eta),
         (T{1}+xi)      * (T{1}-eta) * ((T{0.5625}*p-T{0.0625})*(T{1}-xi+eta) - (T{0.5625}*p+T{0.1875})*xi*eta),
        -(T{1}-eta*eta) * (T{1}+xi)  * ((T{0.5625}*p-T{0.3125})               - (T{0.5625}*p+T{0.1875})*xi),
         (T{1}+xi)      * (T{1}+eta) * ((T{0.5625}*p-T{0.0625})*(T{1}-xi-eta) + (T{0.5625}*p+T{0.1875})*xi*eta),
        -(T{1}-xi*xi)   * (T{1}+eta) * ((T{0.5625}*p-T{0.3125})               - (T{0.5625}*p+T{0.1875})*eta),
         (T{1}-xi)      * (T{1}+eta) * ((T{0.5625}*p-T{0.0625})*(T{1}+xi-eta) - (T{0.5625}*p+T{0.1875})*xi*eta),
        -(T{1}-eta*eta) * (T{1}-xi)  * ((T{0.5625}*p-T{0.3125})               + (T{0.5625}*p+T{0.1875})*xi)
    );
};

}

#endif
#ifndef FINITE_ELEMENT_2D_BASIS_QUARTIC_SERENDIPITY_HPP
#define FINITE_ELEMENT_2D_BASIS_QUARTIC_SERENDIPITY_HPP

// Базисные функции элемента взяты из статьи:
// Астионенко И.А., Козуб Н.А., Литвиненко Е.И., Хомченко А.Н. Управляемые серендиповы поверхности, сохраняющие межэлементную непрерывность.
// В статье был представлен единственный базис, который авторы сочли наиболее адекватным.

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class quartic_serendipity : public geometry_2d<T, rectangle_element_geometry> {
protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;
    
    explicit quartic_serendipity() = default;
    ~quartic_serendipity() override = default;

    // Нумерация узлов на тетричном серендиповом элементе:  12--11--10--9---8
    //                                                      |               |
    //                                                      13              7
    //                                                      |               |
    //                                                      14              6
    //                                                      |               |
    //                                                      15              5
    //                                                      |               |
    //                                                      0---1---2---3---4
    static constexpr std::array<std::array<T, 2>, 16> nodes = { T{-1.0}, T{-1.0},
                                                                T{-0.5}, T{-1.0},
                                                                T{ 0.0}, T{-1.0},
                                                                T{ 0.5}, T{-1.0},
                                                                T{ 1.0}, T{-1.0},
                                                                T{ 1.0}, T{-0.5},
                                                                T{ 1.0}, T{ 0.0},
                                                                T{ 1.0}, T{ 0.5},
                                                                T{ 1.0}, T{ 1.0},
                                                                T{ 0.5}, T{ 1.0},
                                                                T{ 0.0}, T{ 1.0},
                                                                T{-0.5}, T{ 1.0},
                                                                T{-1.0}, T{ 1.0},
                                                                T{-1.0}, T{ 0.5},
                                                                T{-1.0}, T{ 0.0},
                                                                T{-1.0}, T{-0.5} };

    // N_i = 187/1000  (1 + xi_i xi)(1 + eta_i eta)(-500/561 xi^2 + 311/561 xi_i eta_i xi eta - 500/561 eta^2 - 61/561 (xi_i xi + eta_i eta) + 1)(-2 xi_i xi - 2 eta_i eta + 1), xi_i = +-1, eta_i = +-1, i = 0,4,8,12,
    // N_i = -203/800  (1 - xi^2)(1 + eta_i eta)(1 + 4xi_i xi)(-3200/609 xi_i xi - eta_i eta + 1), xi_i = +-1/2, eta_i = +-1/2, i = 1,3,9,11, перестановка xi и eta местами даст базисы для 5,7,13,15
    // N_i = 1141/2000 (1 - xi^2)(1 + eta_i eta)(-4000/1141 xi^2 - 141/1141 eta_i eta + 1), eta_i = +-1, i = 2,10, перестановка xi и eta местами даст базисы для 6 и 14
    static constexpr auto basis = std::make_tuple(
         (T{1}-xi)      * (T{1}-eta)     * (T{ 2.0000}*(xi     +            eta) + T{1.00000}) * (T{561} + T{61}*(xi+eta) - T{500}*(xi*xi + eta*eta) + T{311}*xi*eta) / T{3000},
        -(T{1}-xi*xi)   * (T{1}-eta)     * (T{ 2}/T{3}*xi      + T{0.25375}*eta  + T{0.25375}) * (T{1}   - T{ 2}* xi),
         (T{1}-xi*xi)   * (T{1}-eta)     * (T{-2.0000}*xi*xi   + T{0.07050}*eta  + T{0.57050}),
        -(T{1}-xi*xi)   * (T{1}-eta)     * (T{-2}/T{3}*xi      + T{0.25375}*eta  + T{0.25375}) * (T{1}   + T{ 2}* xi),
         (T{1}+xi)      * (T{1}-eta)     * (T{-2.0000}*(xi     -            eta) + T{1.00000}) * (T{561} - T{61}*(xi-eta) - T{500}*(xi*xi + eta*eta) - T{311}*xi*eta) / T{3000},
        -(T{1}+xi)      * (T{1}-eta*eta) * (T{ 2}/T{3}*eta     - T{0.25375}*xi   + T{0.25375}) * (T{1}   - T{ 2}* eta),
         (T{1}+xi)      * (T{1}-eta*eta) * (T{-2.0000}*eta*eta - T{0.07050}*xi   + T{0.57050}),
        -(T{1}+xi)      * (T{1}-eta*eta) * (T{-2}/T{3}*eta     - T{0.25375}*xi   + T{0.25375}) * (T{1}   + T{ 2}* eta),
         (T{1}+xi)      * (T{1}+eta)     * (T{-2.0000}*(xi     +            eta) + T{1.00000}) * (T{561} - T{61}*(xi+eta) - T{500}*(xi*xi + eta*eta) + T{311}*xi*eta) / T{3000},
        -(T{1}-xi*xi)   * (T{1}+eta)     * (T{-2}/T{3}*xi      - T{0.25375}*eta  + T{0.25375}) * (T{1}   + T{ 2}* xi),
         (T{1}-xi*xi)   * (T{1}+eta)     * (T{-2.0000}*xi*xi   - T{0.07050}*eta  + T{0.57050}),
        -(T{1}-xi*xi)   * (T{1}+eta)     * (T{ 2}/T{3}*xi      - T{0.25375}*eta  + T{0.25375}) * (T{1}   - T{ 2}* xi),
         (T{1}-xi)      * (T{1}+eta)     * (T{ 2.0000}*(xi     -            eta) + T{1.00000}) * (T{561} + T{61}*(xi-eta) - T{500}*(xi*xi + eta*eta) - T{311}*xi*eta) / T{3000},
        -(T{1}-xi)      * (T{1}-eta*eta) * (T{-2}/T{3}*eta     + T{0.25375}*xi   + T{0.25375}) * (T{1}   + T{ 2}* eta),
         (T{1}-xi)      * (T{1}-eta*eta) * (T{-2.0000}*eta*eta + T{0.07050}*xi   + T{0.57050}),
        -(T{1}-xi)      * (T{1}-eta*eta) * (T{ 2}/T{3}*eta     + T{0.25375}*xi   + T{0.25375}) * (T{1}   - T{ 2}* eta)
    );
};

}

#endif
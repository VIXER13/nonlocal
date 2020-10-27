#ifndef FINITE_ELEMENT_QUARTIC_SERENDIPITY_ELEMENT_HPP
#define FINITE_ELEMENT_QUARTIC_SERENDIPITY_ELEMENT_HPP

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
    static constexpr std::array<std::array<T, 2>, 16> nodes = { -1.0, -1.0,
                                                                -0.5, -1.0,
                                                                 0.0, -1.0,
                                                                 0.5, -1.0,
                                                                 1.0, -1.0,
                                                                 1.0, -0.5,
                                                                 1.0,  0.0,
                                                                 1.0,  0.5,
                                                                 1.0,  1.0,
                                                                 0.5,  1.0,
                                                                 0.0,  1.0,
                                                                -0.5,  1.0,
                                                                -1.0,  1.0,
                                                                -1.0,  0.5,
                                                                -1.0,  0.0,
                                                                -1.0, -0.5 };

    // N_i = 187/1000  (1 + xi_i xi)(1 + eta_i eta)(-500/561 xi^2 + 311/561 xi_i eta_i xi eta - 500/561 eta^2 - 61/561 (xi_i xi + eta_i eta) + 1)(-2 xi_i xi - 2 eta_i eta + 1), xi_i = +-1, eta_i = +-1, i = 0,4,8,12,
    // N_i = -203/800  (1 - xi^2)(1 + eta_i eta)(1 + 4xi_i xi)(-3200/609 xi_i xi - eta_i eta + 1), xi_i = +-1/2, eta_i = +-1/2, i = 1,3,9,11, перестановка xi и eta местами даст базисы для 5,7,13,15
    // N_i = 1141/2000 (1 - xi^2)(1 + eta_i eta)(-4000/1141 xi^2 - 141/1141 eta_i eta + 1), eta_i = +-1, i = 2,10, перестановка xi и eta местами даст базисы для 6 и 14
    static constexpr auto basis = std::make_tuple(
         (1.-xi)      * (1.-eta)     * ( 2.000*(xi     +         eta) + 1.00000) * (561. + 61.*(xi+eta) - 500.*(xi*xi + eta*eta) + 311.*xi*eta) / 3000.,
        -(1.-xi*xi)   * (1.-eta)     * ( 2./3.*xi      + 0.25375*eta  + 0.25375) * (1.-2.*xi),
         (1.-xi*xi)   * (1.-eta)     * (-2.000*xi*xi   + 0.0705 *eta  + 0.57050),
        -(1.-xi*xi)   * (1.-eta)     * (-2./3.*xi      + 0.25375*eta  + 0.25375) * (1.+2.*xi),
         (1.+xi)      * (1.-eta)     * (-2.000*(xi     -         eta) + 1.00000) * (561. - 61.*(xi-eta) - 500.*(xi*xi + eta*eta) - 311.*xi*eta) / 3000.,
        -(1.+xi)      * (1.-eta*eta) * ( 2./3.*eta     - 0.25375*xi   + 0.25375) * (1.-2.*eta),
         (1.+xi)      * (1.-eta*eta) * (-2.000*eta*eta - 0.0705 *xi   + 0.5705),
        -(1.+xi)      * (1.-eta*eta) * (-2./3.*eta     - 0.25375*xi   + 0.25375) * (1.+2.*eta),
         (1.+xi)      * (1.+eta)     * (-2.000*(xi     +         eta) + 1.00000) * (561. - 61.*(xi+eta) - 500.*(xi*xi + eta*eta) + 311.*xi*eta) / 3000.,
        -(1.-xi*xi)   * (1.+eta)     * (-2./3.*xi      - 0.25375*eta  + 0.25375) * (1.+2.*xi),
         (1.-xi*xi)   * (1.+eta)     * (-2.000*xi*xi   - 0.0705 *eta  + 0.57050),
        -(1.-xi*xi)   * (1.+eta)     * ( 2./3.*xi      - 0.25375*eta  + 0.25375) * (1.-2.*xi),
         (1.-xi)      * (1.+eta)     * ( 2.000*(xi     -         eta) + 1.00000) * (561. + 61.*(xi-eta) - 500.*(xi*xi + eta*eta) - 311.*xi*eta) / 3000.,
        -(1.-xi)      * (1.-eta*eta) * (-2./3.*eta     + 0.25375*xi   + 0.25375) * (1.+2.*eta),
         (1.-xi)      * (1.-eta*eta) * (-2.000*eta*eta + 0.0705 *xi   + 0.57050),
        -(1.-xi)      * (1.-eta*eta) * ( 2./3.*eta     + 0.25375*xi   + 0.25375) * (1.-2.*eta)
    );
};

}

#endif
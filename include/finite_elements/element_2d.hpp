#ifndef FINITE_ELEMENT_2D_HPP
#define FINITE_ELEMENT_2D_HPP

// Классы стратегий для двумерных элементов.

#include "geometry_2d.hpp"

namespace finite_element {

// Функции формы для треугольных элементов конструировались в барицентрических координатах: L_i, i = 0..2.
// Но для простоты дальнейшей работы и оптимальности вычислений, данные координаты были переведены в декартовы: xi и eta.

template<class Type>
class triangle : public geometry_2d<Type, triangle_element_geometry> {
protected:
    explicit triangle() noexcept = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          2--0
    */
    static inline constexpr std::array<std::array<Type, 2>, 3> nodes = { 1.0, 0.0,
                                                                         0.0, 1.0,
                                                                         0.0, 0.0 };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i, i = 0..2                                                                         
    static inline const std::array<std::function<Type(const Type, const Type)>, 3>
        basic_N    = { [](const Type xi, const Type    ) { return         xi; },
                       [](const Type   , const Type eta) { return        eta; },
                       [](const Type xi, const Type eta) { return 1.0-xi-eta; } },

        basic_Nxi  = { [](const Type   , const Type    ) { return  1.0; },
                       [](const Type   , const Type    ) { return  0.0; },
                       [](const Type   , const Type    ) { return -1.0; } },

        basic_Neta = { [](const Type   , const Type    ) { return  0.0; },
                       [](const Type   , const Type    ) { return  1.0; },
                       [](const Type   , const Type    ) { return -1.0; } };
};

template<class Type>
class quadratic_triangle : public geometry_2d<Type, triangle_element_geometry> {
protected:
    explicit quadratic_triangle() noexcept = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          | \
                                                          4  3
                                                          |    \
                                                          2--5--0
    */
    static inline constexpr std::array<std::array<Type, 2>, 6> nodes = { 1.0, 0.0,
                                                                         0.0, 1.0,
                                                                         0.0, 0.0,
                                                                         0.5, 0.5,
                                                                         0.0, 0.5,
                                                                         0.5, 0.0 };

    // Базисные функции в барицентрических координатах имеют вид: N_i = L_i (2 L_i - 1), i = 0..2,
    //                                                            N_3 = 4 L_0 L_1,
    //                                                            N_4 = 4 L_1 L_2,
    //                                                            N_5 = 4 L_2 L_0.
    static inline const std::array<std::function<Type(const Type, const Type)>, 6>
        basic_N    = { [](const Type xi, const Type    ) { return          xi  * (2.0*xi       - 1.0); },
                       [](const Type   , const Type eta) { return          eta * (2.0*eta      - 1.0); },
                       [](const Type xi, const Type eta) { return (xi+eta-1.0) * (2.0*(xi+eta) - 1.0); },
                       [](const Type xi, const Type eta) { return          4.0 * xi  * eta;            },
                       [](const Type xi, const Type eta) { return         -4.0 * eta * (xi+eta-1.0);   },
                       [](const Type xi, const Type eta) { return         -4.0 * xi  * (xi+eta-1.0);   } },

        basic_Nxi  = { [](const Type xi, const Type    ) { return         4.0 * xi - 1.0; },
                       [](const Type   , const Type    ) { return                    0.0; },
                       [](const Type xi, const Type eta) { return    4.0 * (xi+eta-0.75); },
                       [](const Type   , const Type eta) { return              4.0 * eta; },
                       [](const Type   , const Type eta) { return             -4.0 * eta; },
                       [](const Type xi, const Type eta) { return 4.0 * (1.0-2.0*xi-eta); } },

        basic_Neta = { [](const Type   , const Type    ) { return                    0.0; },
                       [](const Type   , const Type eta) { return        4.0 * eta - 1.0; },
                       [](const Type xi, const Type eta) { return    4.0 * (xi+eta-0.75); },
                       [](const Type xi, const Type    ) { return               4.0 * xi; },
                       [](const Type xi, const Type eta) { return 4.0 * (1.0-xi-2.0*eta); },
                       [](const Type xi, const Type    ) { return              -4.0 * xi; } };
};

template<class Type>
class qubic_triangle : public geometry_2d<Type, triangle_element_geometry> {
protected:
    explicit qubic_triangle() noexcept = default;

    /*
        Нумерация узлов на линейном треугольном элементе: 1\
                                                          5 4
                                                          |  \
                                                          6 9 3
                                                          |    \
                                                          2-7-8-0
    */
    static inline constexpr std::array<std::array<Type, 2>, 10> nodes = {     1.0,     0.0,
                                                                              0.0,     1.0,
                                                                              0.0,     0.0,
                                                                          2.0/3.0, 1.0/3.0,
                                                                          1.0/3.0, 2.0/3.0,
                                                                              0.0, 2.0/3.0,
                                                                              0.0, 1.0/3.0,
                                                                          1.0/3.0,     0.0,
                                                                          2.0/3.0,     0.0,
                                                                          1.0/3.0, 1.0/3.0 };

    // Базисные функции в барицентрических координатах имеют вид: N_i = 0.5 L_i (3 L_i - 1)(3 L_i - 2), i = 0..2,
    //                                                            N_3 = 4.5 L_0 L_1 (3 L_0 - 1),
    //                                                            N_4 = 4.5 L_0 L_1 (3 L_1 - 1),
    //                                                            N_5 = 4.5 L_1 L_2 (3 L_1 - 1),
    //                                                            N_6 = 4.5 L_1 L_2 (3 L_2 - 1),
    //                                                            N_7 = 4.5 L_2 L_0 (3 L_2 - 1),
    //                                                            N_8 = 4.5 L_2 L_0 (3 L_0 - 1),
    //                                                            N_9 =  27 L_0 L_1 L_2
    static inline const std::array<std::function<Type(const Type, const Type)>, 10>
        basic_N    = { [](const Type xi, const Type    ) { return xi           * (1.5*xi      -0.5) * (3.0*xi      -2.0);  },
                       [](const Type   , const Type eta) { return eta          * (1.5*eta     -0.5) * (3.0*eta     -2.0);  },
                       [](const Type xi, const Type eta) { return (1.0-xi-eta) * (1.5*(xi+eta)-0.5) * (3.0*(xi+eta)-2.0);  },
                       [](const Type xi, const Type eta) { return          -xi * eta                * (4.5-13.5*xi);       },
                       [](const Type xi, const Type eta) { return          -xi * eta                * (4.5-13.5*eta);      },
                       [](const Type xi, const Type eta) { return          eta * (xi+eta-1.0)       * (4.5-13.5*eta);      },
                       [](const Type xi, const Type eta) { return         -eta * (xi+eta-1.0)       * (9.0-13.5*(xi+eta)); },
                       [](const Type xi, const Type eta) { return          -xi * (xi+eta-1.0)       * (9.0-13.5*(xi+eta)); },
                       [](const Type xi, const Type eta) { return           xi * (xi+eta-1.0)       * (4.5-13.5*xi);       },
                       [](const Type xi, const Type eta) { return 27.0 *    xi * eta                * (1.0-xi-eta);        } },

        basic_Nxi  = { [](const Type xi, const Type    ) { return                                      (13.5*xi - 9.0)*xi + 1.0; },
                       [](const Type   , const Type    ) { return                                                           0.0; },
                       [](const Type xi, const Type eta) { return     -13.5*(xi*xi+eta*eta) - 27.0*xi*eta + 18.0*(xi+eta) - 5.5; },
                       [](const Type xi, const Type eta) { return                                   ( 27.0*xi       -  4.5)*eta; },
                       [](const Type   , const Type eta) { return                                   ( 13.5*eta      -  4.5)*eta; },
                       [](const Type   , const Type eta) { return                                   (-13.5*eta      +  4.5)*eta; },
                       [](const Type xi, const Type eta) { return                                   ( 27.0*(xi+eta) - 22.5)*eta; },
                       [](const Type xi, const Type eta) { return ( 40.5*xi + 54.0*eta - 45.0)*xi + (13.5*eta - 22.5)*eta + 9.0; },
                       [](const Type xi, const Type eta) { return (-40.5*xi - 27.0*eta + 36.0)*xi +               4.5*eta - 4.5; },
                       [](const Type xi, const Type eta) { return                                (27.0*(1.0-eta) - 54.0*xi)*eta; } },

        basic_Neta = { [](const Type   , const Type    ) { return                                                          0.0; },
                       [](const Type   , const Type eta) { return                                   (13.5*eta - 9.0)*eta + 1.0; },
                       [](const Type xi, const Type eta) { return    -13.5*(xi*xi+eta*eta) - 27.0*xi*eta + 18.0*(xi+eta) - 5.5; },
                       [](const Type xi, const Type    ) { return                                         ( 13.5*xi  - 4.5)*xi; },
                       [](const Type xi, const Type eta) { return                                         ( 27.0*eta - 4.5)*xi; },
                       [](const Type xi, const Type eta) { return (-40.5*eta - 27.0*xi + 36.0)*eta +              4.5*xi - 4.5; },
                       [](const Type xi, const Type eta) { return ( 40.5*eta + 54.0*xi - 45.0)*eta + (13.5*xi - 22.5)*xi + 9.0; },
                       [](const Type xi, const Type eta) { return                                    (27.0*(xi+eta) - 22.5)*xi; },
                       [](const Type xi, const Type    ) { return                                         (-13.5*xi  + 4.5)*xi; },
                       [](const Type xi, const Type eta) { return                                (27.0*(1.0-xi) - 54.0*eta)*xi; } };
};

template<class Type>
class bilinear : public geometry_2d<Type, rectangle_element_geometry> {
protected:
    explicit bilinear() noexcept = default;

    // Нумерация узлов на билинейном элементе: 3---2
    //                                         |   |
    //                                         0---1
    static inline constexpr std::array<std::array<Type, 2>, 4> nodes = { -1.0, -1.0,
                                                                          1.0, -1.0,
                                                                          1.0,  1.0,
                                                                         -1.0,  1.0 };

    // Базисные функции в локальной системе координат имеют вид: N_i = 0.25 (1 + xi_i x)(1 + eta_i eta), xi_i =+-1, eta_i = +-1, i = 0..3
    static inline const std::array<std::function<Type(const Type, const Type)>, 4>
        basic_N    = { [](const Type xi, const Type eta) { return 0.25*(1.0-xi)*(1.0-eta); },
                       [](const Type xi, const Type eta) { return 0.25*(1.0+xi)*(1.0-eta); },
                       [](const Type xi, const Type eta) { return 0.25*(1.0+xi)*(1.0+eta); },
                       [](const Type xi, const Type eta) { return 0.25*(1.0-xi)*(1.0+eta); } },

        basic_Nxi  = { [](const Type   , const Type eta) { return -0.25*(1.0-eta); },
                       [](const Type   , const Type eta) { return  0.25*(1.0-eta); },
                       [](const Type   , const Type eta) { return  0.25*(1.0+eta); },
                       [](const Type   , const Type eta) { return -0.25*(1.0+eta); } },

        basic_Neta = { [](const Type xi, const Type    ) { return -0.25*(1.0-xi); },
                       [](const Type xi, const Type    ) { return -0.25*(1.0+xi); },
                       [](const Type xi, const Type    ) { return  0.25*(1.0+xi); },
                       [](const Type xi, const Type    ) { return  0.25*(1.0-xi); } };
};

template<class Type>
class quadratic_serendip : public geometry_2d<Type, rectangle_element_geometry> {
public:
    // В серендиповой аппроксимации высших порядков возникает проблема с негативизмом стандартного базиса в угловых узлах.
    // Для этого вводится специальный параметр p, который позволяет её избежать.
    // В сущности p является значением интеграла по области элемента от угловой функции. Значение интегралов от промежуточных функций есть 1-p.
    static inline Type p = -1.0 / 3.0; // Значение по умолчанию даёт нам классический вариант квадратичных серендиповых элементов.

protected:
    explicit quadratic_serendip() noexcept = default;

    // Нумерация узлов на квадратичном серендиповом элементе: 6---5---4
    //                                                        |       |
    //                                                        7       3
    //                                                        |       |
    //                                                        0---1---2
    static inline constexpr std::array<std::array<Type, 2>, 8> nodes = { -1.0, -1.0,
                                                                          0.0, -1.0,
                                                                          1.0, -1.0,
                                                                          1.0,  0.0,
                                                                          1.0,  1.0,
                                                                          0.0,  1.0,
                                                                         -1.0,  1.0,
                                                                         -1.0,  0.0 };

    // Базисные функции в локальной системе координат имеют вид:
    // N_i = 0.0625 (1 + xi_i xi)(1 + eta_i eta)[(36p-1)(1 - xi_i xi - eta_i eta) + (36p+3)xi_i xi eta_i eta], xi_i = +-1, eta_i = +-1, i = 0,2,4,6,
    // N_i = 0.0625 (1 -  xi^2)(1 + eta_i eta)[(5-36p) + (36p+3)eta_i eta], eta_i = +-1, i = 1,5,
    // N_i = 0.0625 (1 - eta^2)(1 +  xi_i  xi)[(5-36p) + (36p+3) xi_i  xi],  xi_i = +-1, i = 3,7.
    static inline const std::array<std::function<Type(const Type, const Type)>, 8>
        basic_N    = { [](const Type xi, const Type eta) { return  (1.0-xi)      * (1.0-eta) * ((0.5625*p-0.0625)*(1.0+xi+eta) + (0.5625*p+0.1875)*xi*eta); },
                       [](const Type xi, const Type eta) { return -(1.0-xi*xi)   * (1.0-eta) * ((0.5625*p-0.3125)              + (0.5625*p+0.1875)*eta);    },
                       [](const Type xi, const Type eta) { return  (1.0+xi)      * (1.0-eta) * ((0.5625*p-0.0625)*(1.0-xi+eta) - (0.5625*p+0.1875)*xi*eta); },
                       [](const Type xi, const Type eta) { return -(1.0-eta*eta) * (1.0+xi)  * ((0.5625*p-0.3125)              - (0.5625*p+0.1875)*xi);     },
                       [](const Type xi, const Type eta) { return  (1.0+xi)      * (1.0+eta) * ((0.5625*p-0.0625)*(1.0-xi-eta) + (0.5625*p+0.1875)*xi*eta); },
                       [](const Type xi, const Type eta) { return -(1.0-xi*xi)   * (1.0+eta) * ((0.5625*p-0.3125)              - (0.5625*p+0.1875)*eta);    },
                       [](const Type xi, const Type eta) { return  (1.0-xi)      * (1.0+eta) * ((0.5625*p-0.0625)*(1.0+xi-eta) - (0.5625*p+0.1875)*xi*eta); },
                       [](const Type xi, const Type eta) { return -(1.0-eta*eta) * (1.0-xi)  * ((0.5625*p-0.3125)              + (0.5625*p+0.1875)*xi);     } },

        basic_Nxi  = { [](const Type xi, const Type eta) { return -(1.0-eta)     * (xi * ((1.125*p+0.375)*eta + 1.125*p-0.125) - 0.25*eta); },
                       [](const Type xi, const Type eta) { return  (1.0-eta)     *  xi * ((1.125*p+0.375)*eta + 1.125*p-0.625);             },
                       [](const Type xi, const Type eta) { return -(1.0-eta)     * (xi * ((1.125*p+0.375)*eta + 1.125*p-0.125) + 0.25*eta); },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (xi *  (1.125*p+0.375)     + 0.5);                       },
                       [](const Type xi, const Type eta) { return  (1.0+eta)     * (xi * ((1.125*p+0.375)*eta - 1.125*p+0.125) + 0.25*eta); },
                       [](const Type xi, const Type eta) { return -(1.0+eta)     *  xi * ((1.125*p+0.375)*eta - 1.125*p+0.625);             },
                       [](const Type xi, const Type eta) { return  (1.0+eta)     * (xi * ((1.125*p+0.375)*eta - 1.125*p+0.125) - 0.25*eta); },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (xi *  (1.125*p+0.375)     - 0.5);                       } },

        basic_Neta = { [](const Type xi, const Type eta) { return -(1.0-xi)    * (eta * ((1.125*p+0.375)*xi + 1.125*p-0.125) - 0.25*xi); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (eta *  (1.125*p+0.375)    - 0.5);                      },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (eta * ((1.125*p+0.375)*xi - 1.125*p+0.125) - 0.25*xi); },
                       [](const Type xi, const Type eta) { return -(1.0+xi)    *  eta * ((1.125*p+0.375)*xi - 1.125*p+0.625);            },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (eta * ((1.125*p+0.375)*xi - 1.125*p+0.125) + 0.25*xi); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (eta *  (1.125*p+0.375)    + 0.5);                      },
                       [](const Type xi, const Type eta) { return -(1.0-xi)    * (eta * ((1.125*p+0.375)*xi + 1.125*p-0.125) + 0.25*xi); },
                       [](const Type xi, const Type eta) { return  (1.0-xi)    *  eta * ((1.125*p+0.375)*xi + 1.125*p-0.625);            } };
};

template<class Type>
class qubic_serendip : public geometry_2d<Type, rectangle_element_geometry> {
public:
    // В серендиповой аппроксимации высших порядков возникает проблема с негативизмом стандартного базиса в угловых узлах.
    // Для этого вводится специальный параметр p, который позволяет её избежать.
    // В сущности p является значением интеграла по области элемента от угловой функции. Значение интегралов от промежуточных функций есть (1-4p)/2.
    static inline Type p = -0.5; // Значение по умолчанию даёт нам классический вариант кубических серендиповых элементов.

protected:
    explicit qubic_serendip() noexcept = default;

    // Нумерация узлов на кубическом серендиповом элементе: 9---8---7---6
    //                                                      |           |
    //                                                      10          5
    //                                                      |           |
    //                                                      11          4
    //                                                      |           |
    //                                                      0---1---2---3
    static inline constexpr std::array<std::array<Type, 2>, 12> nodes = {     -1.0,     -1.0,
                                                                          -1.0/3.0,     -1.0,
                                                                           1.0/3.0,     -1.0,
                                                                               1.0,     -1.0,
                                                                               1.0, -1.0/3.0,
                                                                               1.0,  1.0/3.0,
                                                                               1.0,      1.0,
                                                                           1.0/3.0,      1.0,
                                                                          -1.0/3.0,      1.0,
                                                                              -1.0,      1.0,
                                                                              -1.0,  1.0/3.0,
                                                                              -1.0, -1.0/3.0 };

    // Базисные функции в локальной системе координат имеют вид:
    // N_i = 1/32 (1 + xi_i xi)(1 + eta_i eta)[9(xi^2 + eta^2) + (18p+9)(xi_i xi eta_i eta - xi_i xi - eta_i eta) + 18p-1], xi_i = +-1, eta_i = +-1, i = 0,3,6,9,
    // N_i = 9/64 (1 -  xi^2)(1 + eta_i eta)[18  xi_i  xi + (2p+1) eta_i eta - 1 + 2p], xi_i = +-1/3, eta_i = +-1  , i = 1,2,7,8,
    // N_i = 9/64 (1 - eta^2)(1 +  xi_i  xi)[18 eta_i eta + (2p+1)  xi_i  xi - 1 + 2p], xi_i = +-1  , eta_i = +-1/3, i = 4,5,10,11.
    static inline const std::array<std::function<Type(const Type, const Type)>, 12>
        basic_N    = { [](const Type xi, const Type eta) { return  (1.0-xi)    * (1.0-eta)     * (0.28125*(xi*xi+eta*eta + (2.00000*p+1.000000)*(xi*eta+xi+eta)) + 0.56250*p - 0.031250); },
                       [](const Type xi, const Type eta) { return -(1.0-xi*xi) * (1.0-eta)     * (0.84375*xi             + (0.28125*p+0.140625)*eta              + 0.28125*p - 0.140625); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (1.0-eta)     * (0.84375*xi             - (0.28125*p+0.140625)*eta              - 0.28125*p + 0.140625); },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (1.0-eta)     * (0.28125*(xi*xi+eta*eta - (2.00000*p+1.000000)*(xi*eta+xi-eta)) + 0.56250*p - 0.031250); },
                       [](const Type xi, const Type eta) { return -(1.0+xi)    * (1.0-eta*eta) * (0.84375*eta            - (0.28125*p+0.140625)*xi               + 0.28125*p - 0.140625); },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (1.0-eta*eta) * (0.84375*eta            + (0.28125*p+0.140625)*xi               - 0.28125*p + 0.140625); },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (1.0+eta)     * (0.28125*(xi*xi+eta*eta + (2.00000*p+1.000000)*(xi*eta-xi-eta)) + 0.56250*p - 0.031250); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (1.0+eta)     * (0.84375*xi             + (0.28125*p+0.140625)*eta              - 0.28125*p + 0.140625); },
                       [](const Type xi, const Type eta) { return -(1.0-xi*xi) * (1.0+eta)     * (0.84375*xi             - (0.28125*p+0.140625)*eta              + 0.28125*p - 0.140625); },
                       [](const Type xi, const Type eta) { return  (1.0-xi)    * (1.0+eta)     * (0.28125*(xi*xi+eta*eta - (2.00000*p+1.000000)*(xi*eta-xi+eta)) + 0.56250*p - 0.031250); },
                       [](const Type xi, const Type eta) { return  (1.0-xi)    * (1.0-eta*eta) * (0.84375*eta            - (0.28125*p+0.140625)*xi               - 0.28125*p + 0.140625); },
                       [](const Type xi, const Type eta) { return -(1.0-xi)    * (1.0-eta*eta) * (0.84375*eta            + (0.28125*p+0.140625)*xi               + 0.28125*p - 0.140625); } },
                      
        basic_Nxi  = { [](const Type xi, const Type eta) { return -(1.0-eta)     * (((1.1250*p+0.56250)*(eta+1.0) + 0.84375*xi  - 0.56250)*xi + 0.28125*eta*eta - 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-eta)     * (((0.5625*p+0.28125)*eta       + 2.53125*xi  - 0.28125     + 0.56250*p) *xi  - 0.84375); },
                       [](const Type xi, const Type eta) { return  (1.0-eta)     * (((0.5625*p+0.28125)*eta       - 2.53125*xi  - 0.28125     + 0.56250*p) *xi  + 0.84375); },
                       [](const Type xi, const Type eta) { return -(1.0-eta)     * (((1.1250*p+0.56250)*(eta+1.0) - 0.84375*xi  - 0.56250)*xi - 0.28125*eta*eta + 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (( 0.5625*p+0.28125)*xi        - 0.84375*eta + 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (( 0.5625*p+0.28125)*xi        + 0.84375*eta + 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0+eta)     * (((1.1250*p+0.56250)*(eta-1.0) + 0.84375*xi  + 0.56250)*xi + 0.28125*eta*eta - 0.31250); },
                       [](const Type xi, const Type eta) { return -(1.0+eta)     * (((0.5625*p+0.28125)*eta       + 2.53125*xi  + 0.28125     - 0.56250*p) *xi  - 0.84375); },
                       [](const Type xi, const Type eta) { return -(1.0+eta)     * (((0.5625*p+0.28125)*eta       - 2.53125*xi  + 0.28125     - 0.56250*p) *xi  + 0.84375); },
                       [](const Type xi, const Type eta) { return  (1.0+eta)     * (((1.1250*p+0.56250)*(eta-1.0) - 0.84375*xi  + 0.56250)*xi - 0.28125*eta*eta + 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (( 0.5625*p+0.28125)*xi        - 0.84375*eta - 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0-eta*eta) * (( 0.5625*p+0.28125)*xi        + 0.84375*eta - 0.28125);                                 } },

        basic_Neta = { [](const Type xi, const Type eta) { return -(1.0-xi)    * (((1.1250*p+0.56250)*(xi+1.0) + 0.84375*eta - 0.56250)*eta + 0.28125*xi*xi  - 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (( 0.5625*p+0.28125)*eta      + 0.84375*xi  - 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (( 0.5625*p+0.28125)*eta      - 0.84375*xi  - 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (((1.1250*p+0.56250)*(xi-1.0) - 0.84375*eta + 0.56250)*eta - 0.28125*xi*xi  + 0.31250); },
                       [](const Type xi, const Type eta) { return -(1.0+xi)    * (((0.5625*p+0.28125)*xi       - 2.53125*eta + 0.28125      - 0.56250*p)*eta + 0.84375); },
                       [](const Type xi, const Type eta) { return -(1.0+xi)    * (((0.5625*p+0.28125)*xi       + 2.53125*eta + 0.28125      - 0.56250*p)*eta - 0.84375); },
                       [](const Type xi, const Type eta) { return  (1.0+xi)    * (((1.1250*p+0.56250)*(xi-1.0) + 0.84375*eta + 0.56250)*eta + 0.28125*xi*xi  - 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (( 0.5625*p+0.28125)*eta      + 0.84375*xi  + 0.28125);                                 },
                       [](const Type xi, const Type eta) { return  (1.0-xi*xi) * (( 0.5625*p+0.28125)*eta      - 0.84375*xi  + 0.28125);                                 },
                       [](const Type xi, const Type eta) { return -(1.0-xi)    * (((1.1250*p+0.56250)*(xi+1.0) - 0.84375*eta - 0.56250)*eta - 0.28125*xi*xi  + 0.31250); },
                       [](const Type xi, const Type eta) { return  (1.0-xi)    * (((0.5625*p+0.28125)*xi       - 2.53125*eta - 0.28125      + 0.56250*p)*eta + 0.84375); },
                       [](const Type xi, const Type eta) { return  (1.0-xi)    * (((0.5625*p+0.28125)*xi       + 2.53125*eta - 0.28125      + 0.56250*p)*eta - 0.84375); }, };
};

}

#endif
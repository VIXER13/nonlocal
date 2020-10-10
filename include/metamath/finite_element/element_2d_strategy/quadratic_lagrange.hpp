#ifndef FINITE_ELEMENT_QUADRATIC_LAGRANGE_ELEMENT_HPP
#define FINITE_ELEMENT_QUADRATIC_LAGRANGE_ELEMENT_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class quadratic_lagrange : public geometry_2d<T, rectangle_element_geometry> {
protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;

    explicit quadratic_lagrange() = default;
    ~quadratic_lagrange() override = default;

    // Нумерация узлов на квадоратичном лагранжевом элементе: 6---5---4
    //                                                        |       |
    //                                                        7   8   3
    //                                                        |       |
    //                                                        0---1---2
    static constexpr std::array<std::array<T, 2>, 9> nodes = { -1., -1.,
                                                                0., -1.,
                                                                1., -1.,
                                                                1.,  0.,
                                                                1.,  1.,
                                                                0.,  1.,
                                                               -1.,  1.,
                                                               -1.,  0.,
                                                                0.,  0. };

    // Базисные функции в локальной системе координат имеют вид: N_i = 0.25 (1 + xi_i x)(1 + eta_i eta), xi_i =+-1, eta_i = +-1, i = 0..3
    static constexpr auto basis = std::make_tuple(
         0.25 * xi  * eta * (xi    - 1.) * (eta     - 1.),
        -0.50 *       eta * (xi*xi - 1.) * (eta     - 1.),
         0.25 * xi  * eta * (xi    + 1.) * (eta     - 1.),
        -0.50 * xi  *       (xi    + 1.) * (eta*eta - 1.),
         0.25 * xi  * eta * (xi    + 1.) * (eta     + 1.),
        -0.50 *       eta * (xi*xi - 1.) * (eta     + 1.),
         0.25 * xi  * eta * (xi    - 1.) * (eta     + 1.),
        -0.50 * xi  *       (xi    - 1.) * (eta*eta - 1.),
                            (xi*xi - 1.) * (eta*eta - 1.)
    );
};

}

#endif
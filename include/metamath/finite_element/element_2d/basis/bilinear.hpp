#ifndef FINITE_ELEMENT_2D_BASIS_BILINEAR_HPP
#define FINITE_ELEMENT_2D_BASIS_BILINEAR_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class bilinear : public geometry_2d<T, rectangle_element_geometry> {
protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;

    explicit bilinear() = default;
    ~bilinear() override = default;

    // Нумерация узлов на билинейном элементе: 3---2
    //                                         |   |
    //                                         0---1
    static constexpr std::array<std::array<T, 2>, 4> nodes = { T{-1}, T{-1},
                                                               T{ 1}, T{-1},
                                                               T{ 1}, T{ 1},
                                                               T{-1}, T{ 1} };

    // Базисные функции в локальной системе координат имеют вид: N_i = 0.25 (1 + xi_i x)(1 + eta_i eta), xi_i =+-1, eta_i = +-1, i = 0..3
    static constexpr auto basis = std::make_tuple(
        T{0.25} * (T{1} - xi) * (T{1} - eta),
        T{0.25} * (T{1} + xi) * (T{1} - eta),
        T{0.25} * (T{1} + xi) * (T{1} + eta),
        T{0.25} * (T{1} - xi) * (T{1} + eta)
    );
};

}

#endif
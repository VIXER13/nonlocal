#ifndef FINITE_ELEMENT_BILINEAR_ELEMENT_HPP
#define FINITE_ELEMENT_BILINEAR_ELEMENT_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class bilinear : public geometry_2d<T, rectangle_element_geometry> {
protected:
    using geometry_2d<T, rectangle_element_geometry>::xi;
    using geometry_2d<T, rectangle_element_geometry>::eta;

    explicit bilinear() = default;

    // Нумерация узлов на билинейном элементе: 3---2
    //                                         |   |
    //                                         0---1
    static constexpr std::array<std::array<T, 2>, 4> nodes = { -1., -1.,
                                                                1., -1.,
                                                                1.,  1.,
                                                               -1.,  1. };

    // Базисные функции в локальной системе координат имеют вид: N_i = 0.25 (1 + xi_i x)(1 + eta_i eta), xi_i =+-1, eta_i = +-1, i = 0..3
    static constexpr auto basis = std::make_tuple(
        0.25 * (1. - xi) * (1. - eta),
        0.25 * (1. + xi) * (1. - eta),
        0.25 * (1. + xi) * (1. + eta),
        0.25 * (1. - xi) * (1. + eta)
    );

    static inline const std::array<std::function<T(const std::array<T, 2>&)>, 4>
        N    = symdiff::to_function<T, 2>(basis),
        Nxi  = symdiff::to_function<T, 2>(symdiff::derivative<xi>(basis)),
        Neta = symdiff::to_function<T, 2>(symdiff::derivative<eta>(basis));
};

}

#endif
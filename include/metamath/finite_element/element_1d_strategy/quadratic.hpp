#ifndef FINITE_ELEMENT_QUADRATIC_ELEMENT_HPP
#define FINITE_ELEMENT_QUADRATIC_ELEMENT_HPP

#include <functional>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class T>
class quadratic : protected geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit quadratic() noexcept = default;

    // Нумерация узлов на квадратичном элементе: 0--1--2
    static constexpr std::array<T, 3> nodes = { -1., 0., 1. };

    static constexpr auto basis = std::make_tuple(
        0.5 * xi * (xi - 1.),
        1.0 - xi * xi,
        0.5 * xi * (xi + 1.)
    );

    static inline const std::array<std::function<T(const std::array<T, 1>&)>, 3>
        N   = symdiff::to_function<T, 1>(basis),
        Nxi = symdiff::to_function<T, 1>(symdiff::derivative<xi>(basis));
};

}

#endif
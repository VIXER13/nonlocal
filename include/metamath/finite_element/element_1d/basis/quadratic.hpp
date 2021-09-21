#ifndef FINITE_ELEMENT_1D_BASIS_QUADRATIC_ELEMENT_HPP
#define FINITE_ELEMENT_1D_BASIS_QUADRATIC_ELEMENT_HPP

#include "geometry_1d.hpp"
#include <functional>

namespace metamath::finite_element {

template<class T>
class quadratic : protected geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit quadratic() noexcept = default;
    ~quadratic() override = default;

    // Нумерация узлов на квадратичном элементе: 0--1--2
    static constexpr std::array<T, 3> nodes = { T{-1}, T{0}, T{1} };

    static constexpr auto basis = std::make_tuple(
        T{0.5} * xi * (xi - T{1}),
        T{1}   - xi *  xi,
        T{0.5} * xi * (xi + T{1})
    );
};

}

#endif
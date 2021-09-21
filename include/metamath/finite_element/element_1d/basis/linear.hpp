#ifndef FINITE_ELEMENT_1D_BASIS_LINEAR_ELEMENT_HPP
#define FINITE_ELEMENT_1D_BASIS_LINEAR_ELEMENT_HPP

#include "geometry_1d.hpp"
#include <functional>

namespace metamath::finite_element {

template<class T>
class linear : public geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit linear() noexcept = default;
    ~linear() override = default;

    // Нумерация узлов на линейном элементе: 0--1
    static constexpr std::array<T, 2> nodes = { T{-1}, T{1} };

    static constexpr auto basis = std::make_tuple(
        T{0.5} * (T{1} - xi),
        T{0.5} * (T{1} + xi)
    );
};

}

#endif
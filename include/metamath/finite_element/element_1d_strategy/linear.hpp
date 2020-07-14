#ifndef FINITE_ELEMENT_LINEAR_ELEMENT_HPP
#define FINITE_ELEMENT_LINEAR_ELEMENT_HPP

#include <functional>
#include "geometry_1d.hpp"

namespace metamath::finite_element {

template<class T>
class linear : public geometry_1d<T, standart_segment_geometry> {
protected:
    using geometry_1d<T, standart_segment_geometry>::xi;

    explicit linear() noexcept = default;

    // Нумерация узлов на линейном элементе: 0--1
    static constexpr std::array<T, 2> nodes = { -1., 1. };

    static constexpr auto basis = std::make_tuple(
        0.5 * (1. - xi),
        0.5 * (1. + xi)
    );

    static inline const std::array<std::function<T(const std::array<T, 1>&)>, 2>
        N   = symdiff::to_function<T, 1>(basis),
        Nxi = symdiff::to_function<T, 1>(symdiff::derivative<xi>(basis));
};

}

#endif
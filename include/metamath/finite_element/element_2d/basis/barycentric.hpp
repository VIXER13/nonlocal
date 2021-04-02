#ifndef FINITE_ELEMENT_2D_BASIS_BARYCENTRIC_HPP
#define FINITE_ELEMENT_2D_BASIS_BARYCENTRIC_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class T>
class barycentric : public geometry_2d<T, triangle_element_geometry> {
protected:
    using geometry_2d<T, triangle_element_geometry>::xi;
    using geometry_2d<T, triangle_element_geometry>::eta;

    explicit barycentric() = default;
    ~barycentric() override = default;

    static constexpr auto L1 = xi;
    static constexpr auto L2 = eta;
    static constexpr auto L3 = T{1} - xi - eta;
};

}

#endif
#ifndef FINITE_ELEMENT_BARYCENTRIC_HPP
#define FINITE_ELEMENT_BARYCENTRIC_HPP

#include "geometry_2d.hpp"

namespace metamath::finite_element {

template<class Type>
class barycentric : public geometry_2d<Type, triangle_element_geometry> {
protected:
    using geometry_2d<Type, triangle_element_geometry>::xi;
    using geometry_2d<Type, triangle_element_geometry>::eta;

    explicit barycentric() = default;

    static constexpr auto L1 = xi;
    static constexpr auto L2 = eta;
    static constexpr auto L3 = 1. - xi - eta;
};

}

#endif
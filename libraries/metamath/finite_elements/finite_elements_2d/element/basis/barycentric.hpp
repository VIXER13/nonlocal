#pragma once

#include "geometry_2d.hpp"
#include "geometric_primitives/triangle.hpp"

namespace metamath::finite_element {

template<class T>
class barycentric : public geometry_2d<T, triangle_element_geometry> {
protected:
    using geometry_2d<T, triangle_element_geometry>::x;
    using geometry_2d<T, triangle_element_geometry>::y;

    static inline constexpr auto L1 = x;
    static inline constexpr auto L2 = y;
    static inline constexpr auto L3 = symbolic::integral_constant<1>{} - x - y;

    explicit barycentric() = default;
    ~barycentric() override = default;
};

}
#pragma once

#include <metamath/finite_elements/finite_elements_2d/geometry/geometry_2d.hpp>
#include <metamath/finite_elements/finite_elements_2d/geometry/geometric_primitives/triangle.hpp>

namespace metamath::finite_element {

template<std::floating_point T, size_t Order>
class barycentric_quadrature;

template<std::floating_point T>
class barycentric_quadrature<T, 1> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 1> nodes = { T{1} / T{3}, T{1} / T{3} };
    static inline constexpr std::array<T, 1> weights = { T{0.5} };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;  
};

template<std::floating_point T>
class barycentric_quadrature<T, 2> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 3> nodes = { T{1} / T{6}, T{1} / T{6},
                                                                      T{2} / T{3}, T{1} / T{6},
                                                                      T{1} / T{6}, T{2} / T{3} };
    static inline constexpr std::array<T, 3> weights = { T{1} / T{6}, T{1} / T{6}, T{1} / T{6} };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;  
};

template<std::floating_point T>
class barycentric_quadrature<T, 3> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 4> nodes = { T{1} / T{3}, T{1} / T{3},
                                                                           T{0.6},      T{0.2},
                                                                           T{0.2},      T{0.6},
                                                                           T{0.2},      T{0.2} };
    static inline constexpr std::array<T, 4> weights = { -T{27} / T{96}, T{25} / T{96}, T{25} / T{96}, T{25} / T{96} };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;
};

template<std::floating_point T>
class barycentric_quadrature<T, 4> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 6> nodes = { T{0.445948490915965}, T{0.445948490915965},
                                                                      T{0.445948490915965}, T{0.108103018168070},
                                                                      T{0.108103018168070}, T{0.445948490915965},
                                                                      T{0.091576213509771}, T{0.091576213509771},
                                                                      T{0.091576213509771}, T{0.816847572980459},
                                                                      T{0.816847572980459}, T{0.091576213509771} };
    static inline constexpr std::array<T, 6> weights = { T{0.11169079483900575}, T{0.11169079483900575}, T{0.11169079483900575},
                                                         T{0.0549758718276611 }, T{0.0549758718276611 }, T{0.0549758718276611 } };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;
};

template<std::floating_point T>
class barycentric_quadrature<T, 5> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 7> nodes = {          T{1} / T{3},          T{1} / T{3},
                                                                      T{0.059715871789770}, T{0.470142064105115},
                                                                      T{0.470142064105115}, T{0.059715871789770},
                                                                      T{0.470142064105115}, T{0.470142064105115},
                                                                      T{0.797426985353087}, T{0.101286507323456},
                                                                      T{0.101286507323456}, T{0.797426985353087},
                                                                      T{0.101286507323456}, T{0.101286507323456} };
    static inline constexpr std::array<T, 7> weights = { T{0.1125},
                                                         T{0.066197076394253 }, T{0.066197076394253 }, T{0.066197076394253 },
                                                         T{0.0629695902724135}, T{0.0629695902724135}, T{0.0629695902724135} };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;
};

template<std::floating_point T>
class barycentric_quadrature<T, 6> : public geometry_2d<T, triangle_element_geometry> {
protected:
    static inline constexpr std::array<std::array<T, 2>, 12> nodes = { T{0.2492867451709105}, T{0.2492867451709105},
                                                                       T{0.2492867451709105}, T{0.501426509658179 },
                                                                       T{0.501426509658179 }, T{0.2492867451709105},
                                                                       T{0.063089014491502 }, T{0.063089014491502 },
                                                                       T{0.063089014491502 }, T{0.873821971016996 },
                                                                       T{0.873821971016996 }, T{0.063089014491502 },
                                                                       T{0.310352475119575 }, T{0.636502475035608 },
                                                                       T{0.636502475035608 }, T{0.053145049844816 },
                                                                       T{0.053145049844816 }, T{0.310352475119575 },
                                                                       T{0.636502475035608 }, T{0.310352475119575 },
                                                                       T{0.310352475119575 }, T{0.053145049844816 } };
    static inline constexpr std::array<T, 12> weights = { 0.0583931378631895, 0.0583931378631895, 0.0583931378631895,
                                                          0.0254224531851035, 0.0254224531851035, 0.0254224531851035,
                                                          0.041425537809187,  0.041425537809187,  0.041425537809187,
                                                          0.041425537809187,  0.041425537809187,  0.041425537809187 };

    explicit barycentric_quadrature() = default;
    ~barycentric_quadrature() override = default;
};

}
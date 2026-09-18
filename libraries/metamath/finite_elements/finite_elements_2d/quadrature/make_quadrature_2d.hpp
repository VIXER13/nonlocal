#pragma once

#include "quadrature_2d.hpp"
#include "primitives/barycentric.hpp"

#include <metamath/finite_elements/finite_elements_1d/quadrature/make_quadrature_1d.hpp>

namespace metamath::finite_element {

enum class geometry_t : uint8_t {
    Rectangle,
    Triangle
};

template<std::floating_point T, size_t Order>
std::unique_ptr<quadrature_2d_base<T>> make_barycentric_quadrature_2d() {
    return std::make_unique<quadrature_2d<T, barycentric_quadrature, Order>>();
}

template<std::floating_point T>
std::unique_ptr<quadrature_2d_base<T>> make_barycentric_quadrature_2d(const size_t order) {
    switch(order) {
        case 1: return make_barycentric_quadrature_2d<T, 1>();
        case 2: return make_barycentric_quadrature_2d<T, 2>();
        case 3: return make_barycentric_quadrature_2d<T, 3>();
        case 4: return make_barycentric_quadrature_2d<T, 4>();
        case 5: return make_barycentric_quadrature_2d<T, 5>();
        case 6: return make_barycentric_quadrature_2d<T, 6>();
        default: throw std::domain_error{"Unsupported quadrature order " + std::to_string(order)};
    }
}

template<std::floating_point T>
std::unique_ptr<quadrature_2d_base<T>> make_quadrature_2d(const size_t order_x, const size_t order_y) {
    return std::make_unique<quadrature_2d<T, cartesian_production>>(*make_quadrature_1d<T>(order_x), *make_quadrature_1d<T>(order_y));
}

template<std::floating_point T>
std::unique_ptr<quadrature_2d_base<T>> make_quadrature_2d(const geometry_t geometry, const size_t order) {
    if (geometry == geometry_t::Triangle)
        return make_barycentric_quadrature_2d<T>(order);
    else if (geometry == geometry_t::Rectangle) {
        const auto quadrature_1d = make_quadrature_1d<T>(order);
        return std::make_unique<quadrature_2d<T, cartesian_production>>(*quadrature_1d, *quadrature_1d);
    } else
        throw std::domain_error{"Unsupported geometry type"};
}

}
#pragma once

#include "cartesian_production.hpp"

#include <metamath/finite_elements/finite_elements_1d/quadrature/make_quadrature_1d.hpp>

namespace metamath::finite_element {

template<std::floating_point T>
std::unique_ptr<quadrature_2d_base<T>> make_quadrature_2d(const size_t order_x, const size_t order_y) {
    return std::make_unique<quadrature_2d<T, cartesian_production>>(*make_quadrature_1d<T>(order_x), *make_quadrature_1d<T>(order_y));
}

template<std::floating_point T>
std::unique_ptr<quadrature_2d_base<T>> make_quadrature_2d(const size_t order) {
    const auto quadrature_1d = make_quadrature_1d<T>(order);
    return std::make_unique<quadrature_2d<T, cartesian_production>>(*quadrature_1d, *quadrature_1d);
}

}
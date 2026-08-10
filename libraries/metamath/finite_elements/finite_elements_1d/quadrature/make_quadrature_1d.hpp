#pragma once

#include "gaussian_quadrature.hpp"
#include "quadrature_1d.hpp"

namespace metamath::finite_element {

template<std::floating_point T, size_t Order>
std::unique_ptr<quadrature_1d_base<T>> make_quadrature_1d() {
    return std::make_unique<quadrature_1d<T, gauss, Order>>();
}

template<std::floating_point T>
std::unique_ptr<quadrature_1d_base<T>> make_quadrature_1d(const size_t order) {
    switch(order) {
        case 1: return make_quadrature_1d<T, 1>();
        case 2: return make_quadrature_1d<T, 2>();
        case 3: return make_quadrature_1d<T, 3>();
        case 4: return make_quadrature_1d<T, 4>();
        case 5: return make_quadrature_1d<T, 5>();
        default: throw std::domain_error{"Unsupported quadrature order " + std::to_string(order)};
    }
}

}
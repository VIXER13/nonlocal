#pragma once

#include <metamath/finite_elements/finite_elements_1d/quadrature/make_quadrature_1d.hpp>

#include "element_1d_integrate.hpp"
#include "element_1d.hpp"
#include "basis/lagrangian_elements_1d.hpp"

namespace metamath::finite_element {

template<std::floating_point T, size_t Order>
std::unique_ptr<element_1d_base<T>> make_element_1d() {
    return std::make_unique<element_1d<T, lagrangian_element_1d, Order>>();
}

template<std::floating_point T>
std::unique_ptr<element_1d_base<T>> make_element_1d(const size_t order) {
    switch(order) {
        case 1: return make_element_1d<T, 1>();
        case 2: return make_element_1d<T, 2>();
        case 3: return make_element_1d<T, 3>();
        case 4: return make_element_1d<T, 4>();
        case 5: return make_element_1d<T, 5>();
        default: throw std::domain_error{"Unsupported element order " + std::to_string(order)};
    }
}

template<std::floating_point T>
element_1d_integrate<T> make_element_1d_integrated(const size_t element_order, const size_t quadrature_order) {
    return element_1d_integrate<T>{*make_element_1d<T>(element_order), *make_quadrature_1d<T>(quadrature_order)};
}

}
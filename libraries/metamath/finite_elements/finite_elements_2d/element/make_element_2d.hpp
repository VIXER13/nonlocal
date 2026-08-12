#pragma once

#include "element_2d_integrate.hpp"
#include "basis/basis_2d.hpp"

#include <metamath/finite_elements/finite_elements_2d/quadrature/make_quadrature_2d.hpp>

namespace metamath::finite_element {

enum class element_t : uint8_t {
    Triangle,
    Serendipity,
    Lagrangian
};

template<std::floating_point T, size_t Order>
std::unique_ptr<element_2d_base<T>> make_triangle_element_2d() {
    return std::make_unique<element_2d<T, triangle, Order>>();
}

template<std::floating_point T>
std::unique_ptr<element_2d_base<T>> make_triangle_element_2d(const size_t order) {
    switch(order) {
        case 1: return make_triangle_element_2d<T, 1>();
        case 2: return make_triangle_element_2d<T, 2>();
        case 3: return make_triangle_element_2d<T, 3>();
        default: throw std::domain_error{"Unsupported triangle element order " + std::to_string(order)};
    }
}

template<std::floating_point T, size_t Order>
std::unique_ptr<element_2d_base<T>> make_serendipity_element_2d() {
    return std::make_unique<element_2d<T, serendipity, Order>>();
}

template<std::floating_point T>
std::unique_ptr<element_2d_base<T>> make_serendipity_element_2d(const size_t order) {
    switch(order) {
        case 1: return make_serendipity_element_2d<T, 1>();
        // Second and third order shall be created manually, due to a bug in compiling specialized templates
        case 2: return std::make_unique<element_2d<T, serendipity, 2>>();
        case 3: return std::make_unique<element_2d<T, serendipity, 3>>();
        case 4: return make_serendipity_element_2d<T, 4>();
        case 5: return make_serendipity_element_2d<T, 5>();
        default: throw std::domain_error{"Unsupported serendipity element order " + std::to_string(order)};
    }
}

template<std::floating_point T, size_t Order_X, size_t Order_Y>
std::unique_ptr<element_2d_base<T>> make_lagrangian_element_2d() {
    return std::make_unique<element_2d<T, lagrangian_element_2d, Order_X, Order_Y>>();
}

template<std::floating_point T>
std::unique_ptr<element_2d_base<T>> make_lagrangian_element_2d(const size_t order) {
    switch(order) {
        case 1: return make_lagrangian_element_2d<T, 1, 1>();
        case 2: return make_lagrangian_element_2d<T, 2, 2>();
        case 3: return make_lagrangian_element_2d<T, 3, 3>();
        case 4: return make_lagrangian_element_2d<T, 4, 4>();
        case 5: return make_lagrangian_element_2d<T, 5, 5>();
        default: throw std::domain_error{"Unsupported lagrangian element order " + std::to_string(order)};
    }
}

template<std::floating_point T>
std::unique_ptr<element_2d_base<T>> make_element_2d(const element_t type, const size_t order) {
    switch(type) {
        case element_t::Triangle: return make_triangle_element_2d<T>(order);
        case element_t::Serendipity: return make_serendipity_element_2d<T>(order);
        case element_t::Lagrangian: return make_lagrangian_element_2d<T>(order);
        default: throw std::domain_error{"Unsupported element type " + std::to_string(static_cast<uint8_t>(type))};
    }
}

template<std::floating_point T>
element_2d_integrate<T> make_element_2d_integrated(const element_t type, const size_t element_order, const size_t quadrature_order) {
    if (type == element_t::Triangle)
        return element_2d_integrate<T>{*make_element_2d<T>(type, element_order), *make_quadrature_2d<T>(geometry_t::Triangle, quadrature_order)};
    return element_2d_integrate<T>{*make_element_2d<T>(type, element_order), *make_quadrature_2d<T>(geometry_t::Rectangle, quadrature_order)};
}

}
#include "init_elements.hpp"

namespace {

using namespace metamath::finite_element;

template<class T, size_t Element_Order_X, size_t Element_Order_Y, size_t Quadrature_Order_X, size_t Quadrature_Order_Y>
std::unique_ptr<element_2d_integrate_base<T>> make_element() {
    return std::make_unique<element_2d_integrate<T, lagrangian_element_2d, Element_Order_X, Element_Order_Y>>(
        quadrature_1d<T, gauss, Quadrature_Order_X>{},
        quadrature_1d<T, gauss, Quadrature_Order_Y>{}
    );
}

template<class T>
std::vector<std::unique_ptr<element_2d_integrate_base<T>>> init() {
    std::vector<std::unique_ptr<element_2d_integrate_base<T>>> result;
    result.emplace_back(make_element<T, 0, 0, 1, 1>());
    result.emplace_back(make_element<T, 0, 1, 1, 1>());
    result.emplace_back(make_element<T, 1, 0, 1, 1>());
    result.emplace_back(make_element<T, 1, 1, 1, 1>());
    result.emplace_back(make_element<T, 1, 2, 1, 2>());
    result.emplace_back(make_element<T, 2, 1, 2, 1>());
    result.emplace_back(make_element<T, 2, 2, 2, 2>());
    result.emplace_back(make_element<T, 2, 3, 2, 2>());
    result.emplace_back(make_element<T, 3, 2, 2, 2>());
    result.emplace_back(make_element<T, 3, 3, 2, 2>());
    return result;
}

}

namespace unit_tests {

template<>
std::vector<std::unique_ptr<element_2d_integrate_base<double>>> init_lagrangian_elements_2d() {
    return init<double>();
}

}
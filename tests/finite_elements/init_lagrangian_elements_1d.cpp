#include "init_elements.hpp"

namespace {

using namespace metamath::finite_element;

template<class T, size_t Element_Order, size_t Quadrature_Order>
std::unique_ptr<element_1d_integrate_base<T>> make_element() {
    return std::make_unique<element_1d_integrate<T, lagrangian_element_1d, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

template<class T>
std::vector<std::unique_ptr<element_1d_integrate_base<T>>> init() {
    std::vector<std::unique_ptr<element_1d_integrate_base<T>>> result;
    result.emplace_back(make_element<T, 0, 1>());
    result.emplace_back(make_element<T, 1, 1>());
    result.emplace_back(make_element<T, 2, 2>());
    result.emplace_back(make_element<T, 3, 2>());
    result.emplace_back(make_element<T, 4, 3>());
    result.emplace_back(make_element<T, 5, 3>());
    return result;
}

}

namespace unit_tests {

template<>
std::vector<std::unique_ptr<element_1d_integrate_base<double>>> init_lagrangian_elements_1d() {
    return init<double>();
}

}
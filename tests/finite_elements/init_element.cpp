#include "init_element.hpp"

namespace {

using namespace metamath::finite_element;

template<class T, size_t Element_Order, size_t Quadrature_Order>
std::unique_ptr<element_1d_integrate_base<T>> make_element() {
    return std::make_unique<element_1d_integrate<T, lagrangian_element_1d, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

template<class T>
std::array<std::unique_ptr<element_1d_integrate_base<T>>, 6> init() {
    return {
        make_element<T, 0, 1>(),
        make_element<T, 1, 1>(),
        make_element<T, 2, 2>(),
        make_element<T, 3, 2>(),
        make_element<T, 4, 3>(),
        make_element<T, 5, 3>()
    };
}

}

namespace unit_tests {

template<>
std::array<std::unique_ptr<element_1d_integrate_base<float>>, 6> init_elements() {
    return init<float>();
}

template<>
std::array<std::unique_ptr<element_1d_integrate_base<double>>, 6> init_elements() {
    return init<double>();
}

template<>
std::array<std::unique_ptr<element_1d_integrate_base<long double>>, 6> init_elements() {
    return init<long double>();
}

}
#include "init_elements.hpp"

namespace {

using namespace metamath::finite_element;

template<class T, size_t Element_Order, size_t Quadrature_Order>
std::unique_ptr<element_2d_integrate_base<T>> make_element() {
    return std::make_unique<element_2d_integrate<T, serendipity, Element_Order>>(
        quadrature_1d<T, gauss, Quadrature_Order>{},
        quadrature_1d<T, gauss, Quadrature_Order>{}
    );
}

template<class T>
std::vector<std::unique_ptr<element_2d_integrate_base<T>>> init() {
    std::vector<std::unique_ptr<element_2d_integrate_base<T>>> result;
    result.emplace_back(make_element<T, 0, 1>());
    result.emplace_back(make_element<T, 1, 1>());
    result.emplace_back(std::make_unique<element_2d_integrate<T, serendipity, 2>>(quadrature_1d<T, gauss, 2>{}, quadrature_1d<T, gauss, 2>{}));
    result.emplace_back(std::make_unique<element_2d_integrate<T, serendipity, 3>>(quadrature_1d<T, gauss, 2>{}, quadrature_1d<T, gauss, 2>{}));
    result.emplace_back(make_element<T, 4, 3>());
    result.emplace_back(make_element<T, 5, 3>());
    return result;
}

}

namespace unit_tests {

template<>
std::vector<std::unique_ptr<element_2d_integrate_base<double>>> init_serendipity_elements_2d() {
    return init<double>();
}

}
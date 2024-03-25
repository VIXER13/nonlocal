#ifndef NONLOCFEM_MAKE_ELEMENT_1D_HPP
#define NONLOCFEM_MAKE_ELEMENT_1D_HPP

#include "metamath.hpp"

using T = double;
using I = int64_t;
static constexpr T epsilon = T{1e-8};

namespace nonstat_1d_tests {

template<class T>
using quadrature_1d_ptr = std::unique_ptr<metamath::finite_element::quadrature_1d_base<T>>;
template<class T>
using finite_element_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

class _make_element_1d final {
    template<class T, size_t Order>
    using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, Order>;
    template<class T, size_t Order>
    using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, Order>;

    explicit constexpr _make_element_1d() noexcept = default;

public:
    template<class T, std::signed_integral I>
    friend quadrature_1d_ptr<T> make_quadrature(const I order);
    template<class T, std::signed_integral I>
    friend finite_element_1d_ptr<T> make_element(const I order, const quadrature_1d_ptr<T>& quadrature);
};

template<class T, std::signed_integral I>
quadrature_1d_ptr<T> make_quadrature(const I order) {
    switch(order) {
        case 1:
            return std::make_unique<_make_element_1d::quadrature<T, 1>>();
        case 2:
            return std::make_unique<_make_element_1d::quadrature<T, 2>>();
        case 3:
            return std::make_unique<_make_element_1d::quadrature<T, 3>>();
        case 4:
            return std::make_unique<_make_element_1d::quadrature<T, 4>>();
        case 5:
            return std::make_unique<_make_element_1d::quadrature<T, 5>>();
        default:
            throw std::logic_error{"Invalid quadrature order " + std::to_string(int(order))};
    }
}

template<class T, std::signed_integral I>
finite_element_1d_ptr<T> make_element(const I order, const quadrature_1d_ptr<T>& quadrature) {
    switch(order) {
        case 1:
            return std::make_unique<_make_element_1d::element_1d<T, 1>>(*quadrature);
        case 2:
            return std::make_unique<_make_element_1d::element_1d<T, 2>>(*quadrature);
        case 3:
            return std::make_unique<_make_element_1d::element_1d<T, 3>>(*quadrature);
        case 4:
            return std::make_unique<_make_element_1d::element_1d<T, 4>>(*quadrature);
        case 5:
            return std::make_unique<_make_element_1d::element_1d<T, 5>>(*quadrature);
        default:
            throw std::logic_error{"Invalid element order " + std::to_string(int(order))};
    }
}

template<class T, std::signed_integral I>
finite_element_1d_ptr<T> make_element(I element_order, I quadrature_order) {
    return make_element(element_order, make_quadrature<T, I>(quadrature_order));
}

}

#endif
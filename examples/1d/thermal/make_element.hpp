#ifndef NONLOCAL_MAKE_ELEMENT_HPP
#define NONLOCAL_MAKE_ELEMENT_HPP

#include "make_quadrature.hpp"

namespace nonlocal {

enum class element_1d_order_t : uint8_t {
    LINEAR = 1,
    QUADRATIC,
    QUBIC,
    QUARTIC,
    QUINTIC
};

template<class T>
using finite_element_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

template<class T>
finite_element_1d_ptr<T> make_element(const element_1d_order_t element_order, 
                                      const quadrature_1d_order_t quadrature_order = quadrature_1d_order_t::DEFAULT);

class _make_element_1d final {
    template<class T, size_t Order>
    using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, Order>;

    explicit constexpr _make_element_1d() noexcept = default;

public:
    template<class T>
    friend finite_element_1d_ptr<T> make_element(const element_1d_order_t element_order, const quadrature_1d_order_t quadrature_order);
};

template<class T>
finite_element_1d_ptr<T> make_element(const element_1d_order_t element_order, const quadrature_1d_order_t quadrature_order) {
    using enum quadrature_1d_order_t;
    switch(element_order) {
        case element_1d_order_t::LINEAR:
            return std::make_unique<_make_element_1d::element_1d<T, 1>>(*make_quadrature<T>(quadrature_order ? : LINEAR));
        case element_1d_order_t::QUADRATIC:
            return std::make_unique<_make_element_1d::element_1d<T, 2>>(*make_quadrature<T>(quadrature_order ? : QUADRATIC));
        case element_1d_order_t::QUBIC:
            return std::make_unique<_make_element_1d::element_1d<T, 3>>(*make_quadrature<T>(quadrature_order ? : QUBIC));
        case element_1d_order_t::QUARTIC:
            return std::make_unique<_make_element_1d::element_1d<T, 4>>(*make_quadrature<T>(quadrature_order ? : QUARTIC));
        case element_1d_order_t::QUINTIC:
            return std::make_unique<_make_element_1d::element_1d<T, 5>>(*make_quadrature<T>(quadrature_order ? : QUINTIC));
        default:
            throw std::logic_error{"Invalid element order " + std::to_string(int(element_order))};
    }
}

}

#endif
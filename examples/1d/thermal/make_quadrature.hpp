#ifndef NONLOCAL_MAKE_QUADRATURE_HPP
#define NONLOCAL_MAKE_QUADRATURE_HPP

#include "metamath.hpp"
#include <memory>

namespace nonlocal {

enum quadrature_1d_order_t : uint8_t {
    DEFAULT,
    LINEAR,
    QUADRATIC,
    QUBIC,
    QUARTIC,
    QUINTIC
};

template<class T>
using quadrature_1d_ptr = std::unique_ptr<metamath::finite_element::quadrature_1d_base<T>>;

class _make_quadrature final {
    template<class T, size_t Order>
    using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, Order>;

    explicit constexpr _make_quadrature() noexcept = default;

public:
    template<class T>
    friend quadrature_1d_ptr<T> make_quadrature(const quadrature_1d_order_t order);
};

template<class T>
quadrature_1d_ptr<T> make_quadrature(const quadrature_1d_order_t order) {
    switch(order) {
        case quadrature_1d_order_t::LINEAR:
            return std::make_unique<_make_quadrature::quadrature<T, 1>>();
        case quadrature_1d_order_t::QUADRATIC:
            return std::make_unique<_make_quadrature::quadrature<T, 2>>();
        case quadrature_1d_order_t::QUBIC:
            return std::make_unique<_make_quadrature::quadrature<T, 3>>();
        case quadrature_1d_order_t::QUARTIC:
            return std::make_unique<_make_quadrature::quadrature<T, 4>>();
        case quadrature_1d_order_t::QUINTIC:
            return std::make_unique<_make_quadrature::quadrature<T, 5>>();
        default:
            throw std::logic_error{"Invalid quadrature order " + std::to_string(int(order))};
    }
}

}

#endif
#ifndef NONLOCFEM_MAKE_ELEMENT_1D_HPP
#define NONLOCFEM_MAKE_ELEMENT_1D_HPP

#include "nonlocal_config.hpp"
#include "metamath.hpp"

using T = double;
using I = int64_t;
static constexpr T epsilon = T{1e-8};

namespace nonstat_1d_tests {

using namespace nonlocal;
using namespace nonlocal::thermal;


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
    template<class T>
    friend quadrature_1d_ptr<T> make_quadrature(const config::order_t order);
    template<class T>
    friend finite_element_1d_ptr<T> make_element(const config::order_t order, const quadrature_1d_ptr<T>& quadrature);
};

template<class T>
quadrature_1d_ptr<T> make_quadrature(const config::order_t order) {
    switch(order) {
        case config::order_t::LINEAR:
            return std::make_unique<_make_element_1d::quadrature<T, 1>>();
        case config::order_t::QUADRATIC:
            return std::make_unique<_make_element_1d::quadrature<T, 2>>();
        case config::order_t::QUBIC:
            return std::make_unique<_make_element_1d::quadrature<T, 3>>();
        case config::order_t::QUARTIC:
            return std::make_unique<_make_element_1d::quadrature<T, 4>>();
        case config::order_t::QUINTIC:
            return std::make_unique<_make_element_1d::quadrature<T, 5>>();
        default:
            throw std::logic_error{"Invalid quadrature order " + std::to_string(int(order))};
    }
}

template<class T>
finite_element_1d_ptr<T> make_element(const config::order_t order, const quadrature_1d_ptr<T>& quadrature) {
    switch(order) {
        case config::order_t::LINEAR:
            return std::make_unique<_make_element_1d::element_1d<T, 1>>(*quadrature);
        case config::order_t::QUADRATIC:
            return std::make_unique<_make_element_1d::element_1d<T, 2>>(*quadrature);
        case config::order_t::QUBIC:
            return std::make_unique<_make_element_1d::element_1d<T, 3>>(*quadrature);
        case config::order_t::QUARTIC:
            return std::make_unique<_make_element_1d::element_1d<T, 4>>(*quadrature);
        case config::order_t::QUINTIC:
            return std::make_unique<_make_element_1d::element_1d<T, 5>>(*quadrature);
        default:
            throw std::logic_error{"Invalid element order " + std::to_string(int(order))};
    }
}

template<class T>
finite_element_1d_ptr<T> make_element(const config::order_t element_order, const config::order_t quadrature_order) {
    return make_element(element_order, make_quadrature<T>(quadrature_order == config::order_t::UNKNOWN ? element_order : quadrature_order));
}

}

#endif
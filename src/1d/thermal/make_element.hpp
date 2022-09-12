#ifndef NONLOCALMPI_MAKE_ELEMENT_HPP
#define NONLOCALMPI_MAKE_ELEMENT_HPP

#include "metamath.hpp"
#include <memory>

namespace nonlocal {

enum class element_type : uint8_t {
    LINEAR = 1,
    QUADRATIC = 2,
    QUBIC = 3,
    QUARTIC = 4,
    QUINTIC = 5
};

template<class T>
using finite_element_1d_ptr = std::unique_ptr<metamath::finite_element::element_1d_integrate_base<T>>;

template<class T>
finite_element_1d_ptr<T> make_element(const element_type type) {
    switch(type) {
        case element_type::LINEAR: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 1>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 1>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUADRATIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 2>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 2>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUBIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 3>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 3>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUARTIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 4>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 4>;
            return std::make_unique<element_1d>(quadrature{});
        }

        case element_type::QUINTIC: {
            using quadrature = metamath::finite_element::quadrature_1d<T, metamath::finite_element::gauss, 5>;
            using element_1d = metamath::finite_element::element_1d_integrate<T, metamath::finite_element::lagrangian_element_1d, 5>;
            return std::make_unique<element_1d>(quadrature{});
        }

        default:
            throw std::logic_error{"Unknown element type " + std::to_string(int(type))};
    }
}

}

#endif
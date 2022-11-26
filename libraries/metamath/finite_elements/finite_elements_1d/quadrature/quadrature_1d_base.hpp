#ifndef FINITE_ELEMENT_QUADRATURE_BASE_1D_HPP
#define FINITE_ELEMENT_QUADRATURE_BASE_1D_HPP

#include "quadrature_base.hpp"
#include "side_1d.hpp"

#include <array>
#include <memory>

namespace metamath::finite_element {

template<class T>
class quadrature_1d_base : public quadrature_base<T> {
public:
    ~quadrature_1d_base() override = default;
    virtual std::unique_ptr<quadrature_1d_base<T>> clone() const = 0;
    virtual const T node(const size_t i) const = 0;
    virtual T boundary(const side_1d bound) const = 0;
};

}

#endif
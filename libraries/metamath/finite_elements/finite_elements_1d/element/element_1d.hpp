#ifndef FINITE_ELEMENT_1D_HPP
#define FINITE_ELEMENT_1D_HPP

#include "element_1d_base.hpp"
#include "derivative_element_basis_1d.hpp"

#include <functional>
#include <tuple>

#include <iostream>

namespace metamath::finite_element {

template<class T, size_t Derivative_Order, template<class, auto...> class Element_Type, auto... Args>
class element_1d : public virtual element_1d_base<T>,
                   public derivative_element_basis_1d<T, 1, Derivative_Order, Element_Type, Args...> {
    using derivative_base = derivative_element_basis_1d<T, 1, Derivative_Order, Element_Type, Args...>;

public:
    using element_1d_base<T>::nodes;

    ~element_1d() override = default;

    size_t nodes_count() const override {
        return derivative_base::basis_and_derivatives.front().size();
    }

    size_t max_derivative_order() const override {
        return derivative_base::basis_and_derivatives.size();
    }

    T node(const size_t i) const override {
        return derivative_base::nodes[i];
    }

    T N(const size_t d, const size_t i, const T xi) const override {
        return derivative_base::basis_and_derivatives[d][i]({xi});
    }
    
    T boundary(const side_1d bound) const override {
        return derivative_base::boundary(bound);
    }
};

}

#endif
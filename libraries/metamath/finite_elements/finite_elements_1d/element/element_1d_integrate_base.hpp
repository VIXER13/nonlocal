#pragma once

#include "element_1d_base.hpp"

#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>


namespace metamath::finite_element {

template<class T>
class element_1d_integrate_base : public element_integrate_base<T>,
                                  private virtual element_1d_base<T> {
protected:
    std::vector<T> _qNxi;
    std::unique_ptr<quadrature_1d_base<T>> _quadrature = nullptr;

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_integrate_base<T>::nodes;
    using element_integrate_base<T>::qnodes;
    using element_1d_base<T>::boundary;
    using element_1d_base<T>::node;
    using element_1d_base<T>::N;
    using element_1d_base<T>::Nxi;

    ~element_1d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature) = 0;
    const quadrature_1d_base<T>& quadrature() const { return *_quadrature; }

    T qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }
};

}
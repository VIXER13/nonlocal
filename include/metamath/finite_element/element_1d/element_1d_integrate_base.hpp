#ifndef FINITE_ELEMENT_1D_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_1D_INTEGRATE_BASE_HPP

#include "element_integrate_base.hpp"
#include "element_1d_base.hpp"
#include "quadrature.hpp"
#include <memory>

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
    using element_1d_base<T>::boundary;
    using element_1d_base<T>::node;
    using element_1d_base<T>::N;
    using element_1d_base<T>::Nxi;

    ~element_1d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature) = 0;
    const std::unique_ptr<quadrature_1d_base<T>>& quadrature() const noexcept { return _quadrature; }

    T qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }
};

}

#endif
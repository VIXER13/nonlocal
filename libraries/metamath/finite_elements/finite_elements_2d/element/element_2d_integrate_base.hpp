#ifndef FINITE_ELEMENT_2D_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_2D_INTEGRATE_BASE_HPP

#include "finite_element_integrate_base.hpp"
#include "element_2d_base.hpp"
#include "quadrature_1d_base.hpp"

namespace metamath::finite_element {

template<class T>
class element_2d_integrate_base : public element_integrate_base<T>,
                                  public virtual element_2d_base<T> {
protected:
    std::vector<T> _qNxi, _qNeta;

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_integrate_base<T>::nodes;
    using element_integrate_base<T>::qnodes;
    using element_2d_base<T>::boundary;
    using element_2d_base<T>::node;
    using element_2d_base<T>::N;
    using element_2d_base<T>::Nxi;
    using element_2d_base<T>::Neta;

    ~element_2d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) = 0;

    T qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    T qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }
};

}

#endif
#ifndef FINITE_ELEMENT_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_INTEGRATE_BASE_HPP

#include <vector>
#include "element_base.hpp"
#include "quadrature.hpp"

namespace metamath::finite_element {

template<class T>
class element_integrate_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    std::vector<T> _weights, _qN;
    explicit element_integrate_base() noexcept = default;

public:
    size_t qnodes_count() const noexcept { return _weights.size(); }
    size_t  nodes_count() const noexcept { return _qN.size() / qnodes_count(); }

    T weight(const size_t q) const noexcept { return _weights[q]; }

    T qN(const size_t i, const size_t q) const noexcept { return _qN[i*qnodes_count() + q]; }

    virtual ~element_integrate_base() noexcept = default;
};

template<class T>
class element_1d_integrate_base : public element_integrate_base<T>,
                                  private virtual element_1d_base<T> {
protected:
    std::vector<T> _qNxi;

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_1d_base<T>::boundary;
    using element_1d_base<T>::node;
    using element_1d_base<T>::N;
    using element_1d_base<T>::Nxi;
    
    virtual void set_quadrature(const quadrature_base<T>& quadrature) = 0;

    T qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }
};

template<class T>
class element_2d_integrate_base : public element_integrate_base<T>,
                                  public virtual element_2d_base<T> {
protected:
    std::vector<T> _qNxi, _qNeta;

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_2d_base<T>::boundary;
    using element_2d_base<T>::node;
    using element_2d_base<T>::N;
    using element_2d_base<T>::Nxi;
    using element_2d_base<T>::Neta;
    
    virtual void set_quadrature(const quadrature_base<T>& quadrature_xi, const quadrature_base<T>& quadrature_eta) = 0;

    T qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    T qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }
};

}

#endif
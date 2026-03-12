#pragma once

#include "element_2d_base.hpp"

#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>
#include <metamath/types/copyable_ptrs.hpp>

namespace metamath::finite_element {

template<class T>
class element_2d_integrate_base : public element_integrate_base<T> {
protected:
    template<class U>
    using cuptr_t = metamath::types::copyable_uptr<U>;
    std::vector<T> _qNxi, _qNeta;
    cuptr_t<element_2d_base<T>> _element;
    

    void set_element(std::unique_ptr<element_2d_base<T>> element) noexcept {
        _element.reset(std::move(element));
    }

    const element_2d_base<T>& element() const noexcept { return *_element; }

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_integrate_base<T>::nodes;
    using element_integrate_base<T>::qnodes;

    ~element_2d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) = 0;

    const std::array<T, 2>& node(const size_t i) const { return element().node(i); }

    T N   (const size_t i, const std::array<T, 2>& x) const { return element().N(i, x); }
    T Nxi (const size_t i, const std::array<T, 2>& x) const { return element().Nxi(i, x); }
    T Neta(const size_t i, const std::array<T, 2>& x) const { return element().Neta(i, x); }

    T boundary(const side_2d bound, const T x) const { return element().boundary(bound, x); }

    T qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    T qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }
};

}
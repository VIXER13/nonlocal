#pragma once

#include "element_1d_base.hpp"

#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>
#include <metamath/types/copyable_ptrs.hpp>


namespace metamath::finite_element {

template<class T>
class element_1d_integrate_base : public element_integrate_base<T>,
                                  public metamath::types::copyable_ptrs<quadrature_1d_base<T>, element_1d_base<T>> {
protected:
    using holder_t = metamath::types::copyable_ptrs<quadrature_1d_base<T>, element_1d_base<T>>;
    using holder_t::_ptrs;
    std::vector<T> _qNxi;

    void set_element(std::unique_ptr<element_1d_base<T>>&& element) noexcept {
        auto& [quadrature_ptr, element_ptr] = _ptrs;
        element_ptr = std::move(element);
    }

    const element_1d_base<T>& element() const noexcept { return *std::get<1>(_ptrs); }

public:
    using element_integrate_base<T>::nodes_count;
    using element_integrate_base<T>::qnodes_count;
    using element_integrate_base<T>::nodes;
    using element_integrate_base<T>::qnodes;

    ~element_1d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature) = 0;
    const quadrature_1d_base<T>& quadrature() const { return *std::get<0>(_ptrs); }

    T node(const size_t i) const { return element().node(i); }
    T N(const size_t i, const T xi) const { return element().N(i, xi); }
    T Nxi(const size_t i, const T xi) const { return element().Nxi(i, xi); }
    T boundary(const side_1d bound) const { return element().boundary(bound); }

    T qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }
};

}
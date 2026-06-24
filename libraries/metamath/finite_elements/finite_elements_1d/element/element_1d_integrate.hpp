#pragma once

#include "element_1d.hpp"

#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>
#include <metamath/types/copyable_ptrs.hpp>

namespace metamath::finite_element {

template<class T>
class element_1d_integrate : public element_integrate_base<T> {
    using element_integrate_base_t = element_integrate_base<T>;
    using element_integrate_base_t::_nearest_qnode;
    using element_integrate_base_t::_weights;
    using element_integrate_base_t::_qN;

    template<class U>
    using cuptr_t = metamath::types::copyable_uptr<U>;

    std::vector<T> _qNxi;
    cuptr_t<quadrature_1d_base<T>> _quadrature;
    cuptr_t<element_1d_base<T>> _element;

    void set_element(std::unique_ptr<element_1d_base<T>> element) noexcept {
        _element.reset(std::move(element));
    }

    const element_1d_base<T>& element() const noexcept { return *_element; }

    element_1d_integrate() = default;

public:
    using element_integrate_base_t::qnodes_count;
    using element_integrate_base_t::nodes_count;
    using element_integrate_base_t::qnodes;
    using element_integrate_base_t::nodes;

    element_1d_integrate(std::unique_ptr<element_1d_base<T>> element, const quadrature_1d_base<T>& quadrature)
    {
        set_element(std::move(element));
        set_quadrature(quadrature);
    }

    ~element_1d_integrate() override = default;

    const quadrature_1d_base<T>& quadrature() const { return *_quadrature; }

    T node(const size_t i) const { return element().node(i); }
    T N(const size_t i, const T xi) const { return element().N(i, xi); }
    T Nxi(const size_t i, const T xi) const { return element().Nxi(i, xi); }
    T boundary(const side_1d bound) const { return element().boundary(bound); }

    T qNxi(const size_t i, const size_t q) const noexcept { return _qNxi[i*qnodes_count() + q]; }

    std::unique_ptr<element_integrate_base<T>> copy() const override {
        return std::make_unique<element_1d_integrate>(*this);
    }

    void set_quadrature(const quadrature_1d_base<T>& quadrature) {
        _quadrature.reset(quadrature.copy());
        std::vector<T> xi(quadrature.nodes_count());
        _weights.resize(quadrature.nodes_count());
        T jacobian = (           boundary(side_1d::RIGHT) -            boundary(side_1d::LEFT)) /
                     (quadrature.boundary(side_1d::RIGHT) - quadrature.boundary(side_1d::LEFT));
        for(size_t q = 0; q < quadrature.nodes_count(); ++q) {
            xi[q] = boundary(side_1d::LEFT) + (quadrature.node(q) - quadrature.boundary(side_1d::LEFT)) * jacobian;
            _weights[q] = quadrature.weight(q) * jacobian;
        }

        _nearest_qnode.resize(element().nodes_count());
        _qN  .resize(element().nodes_count() * qnodes_count());
        _qNxi.resize(element().nodes_count() * qnodes_count());
        for(size_t i = 0; i < nodes_count(); ++i) {
            size_t nearest_quadrature = 0;
            T length = std::numeric_limits<T>::max();
            for(size_t q = 0; q < qnodes_count(); ++q) {
                _qN  [i*qnodes_count() + q] = N  (i, xi[q]);
                _qNxi[i*qnodes_count() + q] = Nxi(i, xi[q]);
                if (const T curr_length = std::abs(xi[q] - node(i)); length > curr_length) {
                    length = curr_length;
                    nearest_quadrature = q;
                }
            }
            _nearest_qnode[i] = nearest_quadrature;
        }
    }
};

}
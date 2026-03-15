#pragma once

#include "element_2d_base.hpp"
#include "element_2d_serendipity.hpp"

#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>
#include <metamath/types/copyable_ptrs.hpp>

namespace metamath::finite_element {

template<class T>
class element_2d_integrate : public element_integrate_base<T> {
protected:
    using element_integrate_base_t = element_integrate_base<T>;
    using element_integrate_base_t::_nearest_qnode;
    using element_integrate_base_t::_weights;
    using element_integrate_base_t::_qN;

    template<class U>
    using cuptr_t = metamath::types::copyable_uptr<U>;

    std::vector<T> _qNxi, _qNeta;
    cuptr_t<element_2d_base<T>> _element;

    void set_element(std::unique_ptr<element_2d_base<T>> element) noexcept {
        _element.reset(std::move(element));
    }

    const element_2d_base<T>& element() const noexcept { return *_element; }

    element_2d_integrate() = default;

public:
    using element_integrate_base_t::qnodes_count;
    using element_integrate_base_t::nodes_count;
    using element_integrate_base_t::nodes;
    using element_integrate_base_t::qnodes;

    const std::array<T, 2>& node(const size_t i) const { return element().node(i); }

    T N   (const size_t i, const std::array<T, 2>& x) const { return element().N(i, x); }
    T Nxi (const size_t i, const std::array<T, 2>& x) const { return element().Nxi(i, x); }
    T Neta(const size_t i, const std::array<T, 2>& x) const { return element().Neta(i, x); }

    T boundary(const side_2d bound, const T x) const { return element().boundary(bound, x); }

    T qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    T qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }

    explicit element_2d_integrate(std::unique_ptr<element_2d_base<T>> element, const quadrature_1d_base<T>& quadrature) {
        set_element(std::move(element));
        set_quadrature(quadrature, quadrature);
    }

    explicit element_2d_integrate(std::unique_ptr<element_2d_base<T>> element,
                                 const quadrature_1d_base<T>& quadrature_x,
                                 const quadrature_1d_base<T>& quadrature_y) {
        set_element(std::move(element));
        set_quadrature(quadrature_x, quadrature_y);
    }

    ~element_2d_integrate() override = default;

    std::unique_ptr<element_integrate_base<T>> copy() const override {
        return std::make_unique<element_2d_integrate>(*this);
    }

    void set_quadrature(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) {
        T jacobian_x = (             boundary(side_2d::RIGHT, 0) -              boundary(side_2d::LEFT, 0)) /
                       (quadrature_x.boundary(side_1d::RIGHT   ) - quadrature_x.boundary(side_1d::LEFT   ));
        std::vector<T> jacobian_y(quadrature_x.nodes_count());
        std::vector<T> x(quadrature_x.nodes_count());
        std::vector<T> y(quadrature_y.nodes_count());

        _weights.resize(quadrature_x.nodes_count() * quadrature_y.nodes_count());
        for(size_t i = 0; i < quadrature_x.nodes_count(); ++i) {
            x[i] = boundary(side_2d::LEFT, 0) + (quadrature_x.node(i) - quadrature_x.boundary(side_1d::LEFT)) * jacobian_x;
            jacobian_y[i] = (             boundary(side_2d::UP,  x[i]) -              boundary(side_2d::DOWN, x[i])) /
                            (quadrature_y.boundary(side_1d::RIGHT    ) - quadrature_y.boundary(side_1d::LEFT      ));
            for(size_t j = 0; j < quadrature_y.nodes_count(); ++j) {
                y[j] = boundary(side_2d::DOWN, x[i]) + (quadrature_y.node(j)-quadrature_y.boundary(side_1d::LEFT)) * jacobian_y[i];
                _weights[i*quadrature_y.nodes_count() + j] = quadrature_x.weight(i) * jacobian_x * quadrature_y.weight(j) * jacobian_y[i];
            }
        }

        _nearest_qnode.resize(element().nodes_count());
        _qN.resize(element().nodes_count() * qnodes_count());
        _qNxi.resize(element().nodes_count() * qnodes_count());
        _qNeta.resize(element().nodes_count() * qnodes_count());
        for(size_t i = 0; i < nodes_count(); ++i) {
            size_t nearest_quadrature = 0;
            T length = std::numeric_limits<T>::max();
            for(size_t j = 0; j < quadrature_x.nodes_count(); ++j)
                for(size_t k = 0; k < quadrature_y.nodes_count(); ++k) {
                    _qN   [i*qnodes_count() + j*quadrature_y.nodes_count() + k] = N   (i, {x[j], y[k]});
                    _qNxi [i*qnodes_count() + j*quadrature_y.nodes_count() + k] = Nxi (i, {x[j], y[k]});
                    _qNeta[i*qnodes_count() + j*quadrature_y.nodes_count() + k] = Neta(i, {x[j], y[k]});
                    if (const T curr_length = functions::distance(node(i), {x[j], y[k]}); length > curr_length) {
                        length = curr_length;
                        nearest_quadrature = j*quadrature_y.nodes_count() + k;
                    }
                }
            _nearest_qnode[i] = nearest_quadrature;
        }
    }
};

}
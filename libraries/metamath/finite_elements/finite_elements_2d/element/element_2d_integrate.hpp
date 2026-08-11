#pragma once

#include "element_2d_serendipity.hpp"

#include <metamath/linear/distance.hpp>
#include <metamath/finite_elements/base/finite_element_integrate_base.hpp>
#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d_base.hpp>
#include <metamath/finite_elements/finite_elements_2d/quadrature/quadrature_2d_base.hpp>
#include <metamath/finite_elements/finite_elements_2d/geometry/geometric_primitives/triangle.hpp>
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
    cuptr_t<quadrature_2d_base<T>> _quadrature;

    element_2d_integrate() = default;

public:
    using element_integrate_base_t::qnodes_count;
    using element_integrate_base_t::nodes_count;
    using element_integrate_base_t::nodes;
    using element_integrate_base_t::qnodes;

    explicit element_2d_integrate(std::unique_ptr<element_2d_base<T>> element, const quadrature_1d_base<T>& quadrature)
        : _element{element->copy()} {
        set_quadrature(quadrature, quadrature);
    }

    explicit element_2d_integrate(std::unique_ptr<element_2d_base<T>> element,
                                 const quadrature_1d_base<T>& quadrature_x,
                                 const quadrature_1d_base<T>& quadrature_y)
        : _element{element->copy()} {
        set_quadrature(quadrature_x, quadrature_y);
    }

    explicit element_2d_integrate(element_2d_base<T>& element, const quadrature_2d_base<T>& quadrature)
        : _element{element.copy()} {
        set_quadrature(quadrature);
    }

    ~element_2d_integrate() override = default;

    const element_2d_base<T>& element() const noexcept { return *_element; }
    const quadrature_2d_base<T>& quadrature() const noexcept { return *_quadrature; }

    T qNxi (const size_t i, const size_t q) const noexcept { return _qNxi [i*qnodes_count() + q]; }
    T qNeta(const size_t i, const size_t q) const noexcept { return _qNeta[i*qnodes_count() + q]; }

    std::unique_ptr<element_integrate_base<T>> copy() const override {
        return std::make_unique<element_2d_integrate>(*this);
    }

    void set_quadrature(const quadrature_2d_base<T>& quadrature) {
        const bool is_triangle_element = bool(dynamic_cast<const triangle_element_geometry<T>*>(&element()));
        const bool is_triangle_quadrature = bool(dynamic_cast<const triangle_element_geometry<T>*>(&quadrature));
        if (is_triangle_element != is_triangle_quadrature) {
            // TODO: support integration triangle elements by rectangle quadrature if it's possible
            throw std::invalid_argument("element and quadrature shall be of the same type");
        }

        _quadrature = quadrature.copy();
        _weights.resize(quadrature().nodes_count());
        for(const size_t i : std::ranges::iota_view(0zu, quadrature().nodes_count()))
            _weights[i] = quadrature().weight(i);

        _nearest_qnode.resize(element().nodes_count(), 0);
        _qN.resize(element().nodes_count() * qnodes_count());
        _qNxi.resize(element().nodes_count() * qnodes_count());
        _qNeta.resize(element().nodes_count() * qnodes_count());
        for(const size_t i : std::ranges::iota_view(0zu, element().nodes_count())) {
            T length = std::numeric_limits<T>::max();
            for(const size_t j : std::ranges::iota_view(0zu, quadrature().nodes_count())) {
                const auto& qnode = quadrature().node(j);
                _qN   [i * qnodes_count() + j] = element().N   (i, qnode);
                _qNxi [i * qnodes_count() + j] = element().Nxi (i, qnode);
                _qNeta[i * qnodes_count() + j] = element().Neta(i, qnode);
                if (const T curr_length = linear::distance(element().node(i), qnode); length > curr_length) {
                    length = curr_length;
                    _nearest_qnode[i] = j;
                }
            }
        }
    } 

    void set_quadrature(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) {
        T jacobian_x = (element().boundary(side_2d::RIGHT, 0) - element().boundary(side_2d::LEFT, 0)) /
                       (quadrature_x.boundary(side_1d::RIGHT) - quadrature_x.boundary(side_1d::LEFT));
        std::vector<T> jacobian_y(quadrature_x.nodes_count());
        std::vector<T> x(quadrature_x.nodes_count());
        std::vector<T> y(quadrature_y.nodes_count());

        _weights.resize(quadrature_x.nodes_count() * quadrature_y.nodes_count());
        for(size_t i = 0; i < quadrature_x.nodes_count(); ++i) {
            x[i] = element().boundary(side_2d::LEFT, 0) + (quadrature_x.node(i) - quadrature_x.boundary(side_1d::LEFT)) * jacobian_x;
            jacobian_y[i] = (element().boundary(side_2d::UP, x[i]) - element().boundary(side_2d::DOWN, x[i])) /
                            (quadrature_y.boundary(side_1d::RIGHT) - quadrature_y.boundary(side_1d::LEFT));
            for(size_t j = 0; j < quadrature_y.nodes_count(); ++j) {
                y[j] = element().boundary(side_2d::DOWN, x[i]) + (quadrature_y.node(j)-quadrature_y.boundary(side_1d::LEFT)) * jacobian_y[i];
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
                    _qN   [i*qnodes_count() + j*quadrature_y.nodes_count() + k] = element().N   (i, {x[j], y[k]});
                    _qNxi [i*qnodes_count() + j*quadrature_y.nodes_count() + k] = element().Nxi (i, {x[j], y[k]});
                    _qNeta[i*qnodes_count() + j*quadrature_y.nodes_count() + k] = element().Neta(i, {x[j], y[k]});
                    if (const T curr_length = linear::distance(element().node(i), {x[j], y[k]}); length > curr_length) {
                        length = curr_length;
                        nearest_quadrature = j*quadrature_y.nodes_count() + k;
                    }
                }
            _nearest_qnode[i] = nearest_quadrature;
        }
    }
};

}
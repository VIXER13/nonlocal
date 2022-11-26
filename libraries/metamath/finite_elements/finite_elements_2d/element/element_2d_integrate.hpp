#ifndef FINITE_ELEMENT_2D_INTEGRATE_HPP
#define FINITE_ELEMENT_2D_INTEGRATE_HPP

#include "element_2d_integrate_base.hpp"
#include "element_2d.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Element_Type, auto... Args>
class element_2d_integrate : public element_2d_integrate_base<T>,
                             public element_2d<T, Element_Type, Args...> {
protected:
    using element_2d_t = element_2d<T, Element_Type, Args...>;
    using element_integrate_2d_t = element_2d_integrate_base<T>;
    using element_integrate_2d_t::_nearest_qnode;
    using element_integrate_2d_t::_weights;
    using element_integrate_2d_t::_qN;
    using element_integrate_2d_t::_qNxi;
    using element_integrate_2d_t::_qNeta;

public:
    using element_integrate_2d_t::qnodes_count;
    using element_integrate_2d_t::nodes_count;
    using element_2d_t::N;
    using element_2d_t::Nxi;
    using element_2d_t::Neta;
    using element_2d_t::node;
    using element_2d_t::boundary;

    explicit element_2d_integrate(const quadrature_1d_base<T>& quadrature) {
        set_quadrature(quadrature, quadrature);
    }

    explicit element_2d_integrate(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) {
        set_quadrature(quadrature_x, quadrature_y);
    }

    ~element_2d_integrate() override = default;

    void set_quadrature(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) override {
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

        _nearest_qnode.resize(element_2d_t::nodes_count());
        _qN.resize(element_2d_t::nodes_count() * qnodes_count());
        _qNxi.resize(element_2d_t::nodes_count() * qnodes_count());
        _qNeta.resize(element_2d_t::nodes_count() * qnodes_count());
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

#endif
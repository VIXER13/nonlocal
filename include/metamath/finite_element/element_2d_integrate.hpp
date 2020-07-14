#ifndef FINITE_ELEMENT_2D_INTEGRATE_HPP
#define FINITE_ELEMENT_2D_INTEGRATE_HPP

#include "element_integrate_base.hpp"
#include "element_2d.hpp"

namespace metamath::finite_element {

template<class T, template<class> class Element_Type>
class element_2d_integrate : public element_2d_integrate_base<T>,
                             public element_2d<T, Element_Type> {
    using element_2d_integrate_base<T>::_weights;
    using element_2d_integrate_base<T>::_qN;
    using element_2d_integrate_base<T>::_qNxi;
    using element_2d_integrate_base<T>::_qNeta;

public:
    using element_2d_integrate_base<T>::qnodes_count;
    using element_2d_integrate_base<T>::nodes_count;
    using element_2d<T, Element_Type>::N;
    using element_2d<T, Element_Type>::Nxi;
    using element_2d<T, Element_Type>::Neta;
    using element_2d<T, Element_Type>::boundary;

    explicit element_2d_integrate(const quadrature_base<T>& quadrature) {
        set_quadrature(quadrature, quadrature);
    }

    explicit element_2d_integrate(const quadrature_base<T>& quadrature_xi, const quadrature_base<T>& quadrature_eta) {
        set_quadrature(quadrature_xi, quadrature_eta);
    }

    // Узлы в которых происходит расчёт получаются путём декартова произведения квадратур с учётом того,
    // что геометрия описывается таким образом, что параметрическую зависимость имеет верхняя и нижняя границы,
    // а правая и левая заданы константами.
    void set_quadrature(const quadrature_base<T>& quadrature_xi, const quadrature_base<T>& quadrature_eta) override
    {
        T jacobian_xi = (              boundary(side_2d::RIGHT, 0) -               boundary(side_2d::LEFT, 0)) /
                        (quadrature_xi.boundary(side_1d::RIGHT   ) - quadrature_xi.boundary(side_1d::LEFT   ));
        std::vector<T> jacobian_eta(quadrature_xi.nodes_count()),
                       xi(quadrature_xi.nodes_count()),
                       eta(quadrature_eta.nodes_count());

        _weights.resize(quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        for(size_t i = 0; i < quadrature_xi.nodes_count(); ++i) {
            xi[i] = boundary(side_2d::LEFT, 0) + (quadrature_xi.node(i) - quadrature_xi.boundary(side_1d::LEFT)) * jacobian_xi;
            jacobian_eta[i] = (               boundary(side_2d::UP, xi[i]) -                boundary(side_2d::DOWN, xi[i])) /
                              (quadrature_eta.boundary(side_1d::RIGHT    ) - quadrature_eta.boundary(side_1d::LEFT       ));
            for(size_t j = 0; j < quadrature_eta.nodes_count(); ++j) {
                eta[j] = boundary(side_2d::DOWN, xi[i]) + (quadrature_eta.node(j)-quadrature_eta.boundary(side_1d::LEFT)) * jacobian_eta[i];
                _weights[i*quadrature_eta.nodes_count() + j] = quadrature_xi.weight(i) * jacobian_xi * quadrature_eta.weight(j) * jacobian_eta[i];
            }
        }

        _qN.resize(element_2d<T, Element_Type>::nodes_count() * qnodes_count());
        _qNxi.resize(element_2d<T, Element_Type>::nodes_count() * qnodes_count());
        _qNeta.resize(element_2d<T, Element_Type>::nodes_count() * qnodes_count());
        for(size_t i = 0; i < nodes_count(); ++i)
            for(size_t j = 0; j < quadrature_xi.nodes_count(); ++j)
                for(size_t k = 0; k < quadrature_eta.nodes_count(); ++k) {
                    _qN   [i*qnodes_count() + j*quadrature_eta.nodes_count() + k] = N   (i, {xi[j], eta[k]});
                    _qNxi [i*qnodes_count() + j*quadrature_eta.nodes_count() + k] = Nxi (i, {xi[j], eta[k]});
                    _qNeta[i*qnodes_count() + j*quadrature_eta.nodes_count() + k] = Neta(i, {xi[j], eta[k]});
                }
    }
};

}

#endif
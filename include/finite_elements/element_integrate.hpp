#ifndef FINITE_ELEMENT_INTEGRATE_HPP
#define FINITE_ELEMENT_INTEGRATE_HPP

#include "element_base.hpp"

namespace finite_element {

template<class Type>
class element_integrate_base {
    static_assert(std::is_floating_point_v<Type>, "The Type must be floating point.");

protected:
    std::vector<Type> weights;
    matrix<Type> N_in_quad_nodes;
    explicit element_integrate_base() noexcept = default;

public:
    size_t qnodes_count() const noexcept { return N_in_quad_nodes.cols(); }
    size_t  nodes_count() const noexcept { return N_in_quad_nodes.rows(); }

    Type weight(const size_t q) const noexcept { return weights[q]; }

    Type qN(const size_t i, const size_t q) const noexcept { return N_in_quad_nodes(i, q); }

    virtual ~element_integrate_base() noexcept = default;
};

template<class Type>
class element_1d_integrate_base : public element_integrate_base<Type>,
                                  private virtual element_1d_base<Type> {
protected:
    matrix<Type> Nxi_in_quad_nodes;

public:
    using element_integrate_base<Type>::nodes_count;
    using element_integrate_base<Type>::qnodes_count;
    using element_1d_base<Type>::boundary;
    using element_1d_base<Type>::node;
    using element_1d_base<Type>::N;
    using element_1d_base<Type>::Nxi;
    
    virtual void set_quadrature(const quadrature_base<Type>& quadrature) = 0;

    Type qNxi(const size_t i, const size_t q) const noexcept { return Nxi_in_quad_nodes(i, q); }
};

template<class Type, template<class> class Element_Type>
class element_1d_integrate : public element_1d_integrate_base<Type>,
                             public element_1d<Type, Element_Type> {
    using element_1d_integrate_base<Type>::weights;
    using element_1d_integrate_base<Type>::N_in_quad_nodes;
    using element_1d_integrate_base<Type>::Nxi_in_quad_nodes;

public:
    using element_1d_integrate_base<Type>::qnodes_count;
    using element_1d_integrate_base<Type>::nodes_count;
    using element_1d<Type, Element_Type>::boundary;
    using element_1d<Type, Element_Type>::N;
    using element_1d<Type, Element_Type>::Nxi;

    explicit element_1d_integrate(const quadrature_base<Type>& quadrature) { set_quadrature(quadrature); }

    void set_quadrature(const quadrature_base<Type>& quadrature) override {
        std::vector<Type> xi(quadrature.nodes_count());
        weights.resize(quadrature.nodes_count());
        Type jacobian = (           boundary(side_1d::RIGHT) -            boundary(side_1d::LEFT)) /
                        (quadrature.boundary(side_1d::RIGHT) - quadrature.boundary(side_1d::LEFT));
        for(size_t q = 0; q < quadrature.nodes_count(); ++q) {
            xi[q] = boundary(side_1d::LEFT) + (quadrature.node(q) - quadrature.boundary(side_1d::LEFT)) * jacobian;
            weights[q] = quadrature.weight(q) * jacobian;
        }

        N_in_quad_nodes  .resize(element_1d<Type, Element_Type>::nodes_count(), quadrature.nodes_count());
        Nxi_in_quad_nodes.resize(element_1d<Type, Element_Type>::nodes_count(), quadrature.nodes_count());
        for(size_t i = 0; i < nodes_count(); ++i) {
            for(size_t q = 0; q < qnodes_count(); ++q) {
                N_in_quad_nodes  (i, q) = N  (i, xi[q]);
                Nxi_in_quad_nodes(i, q) = Nxi(i, xi[q]);
            }
        }
    }
};

template<class Type>
class element_2d_integrate_base : public element_integrate_base<Type>,
                                  public virtual element_2d_base<Type> {
protected:
    matrix<Type> Nxi_in_quad_nodes, Neta_in_quad_nodes;

public:
    using element_integrate_base<Type>::nodes_count;
    using element_integrate_base<Type>::qnodes_count;
    using element_2d_base<Type>::boundary;
    using element_2d_base<Type>::node;
    using element_2d_base<Type>::N;
    using element_2d_base<Type>::Nxi;
    using element_2d_base<Type>::Neta;
    
    virtual void set_quadrature(const quadrature_base<Type>& quadrature_xi, const quadrature_base<Type>& quadrature_eta) = 0;

    Type qNxi (const size_t i, const size_t q) const noexcept { return Nxi_in_quad_nodes (i, q); }
    Type qNeta(const size_t i, const size_t q) const noexcept { return Neta_in_quad_nodes(i, q); }
};

template<class Type, template<class> class Element_Type>
class element_2d_integrate : public element_2d_integrate_base<Type>,
                             public element_2d<Type, Element_Type> {
    using element_2d_integrate_base<Type>::weights;
    using element_2d_integrate_base<Type>::N_in_quad_nodes;
    using element_2d_integrate_base<Type>::Nxi_in_quad_nodes;
    using element_2d_integrate_base<Type>::Neta_in_quad_nodes;

public:
    using element_2d_integrate_base<Type>::qnodes_count;
    using element_2d_integrate_base<Type>::nodes_count;
    using element_2d<Type, Element_Type>::N;
    using element_2d<Type, Element_Type>::Nxi;
    using element_2d<Type, Element_Type>::Neta;
    using element_2d<Type, Element_Type>::boundary;

    explicit element_2d_integrate(const quadrature_base<Type>& quadrature_xi, const quadrature_base<Type>& quadrature_eta) {
        set_quadrature(quadrature_xi, quadrature_eta);
    }

    // Узлы в которых происходит расчёт получаются путём декартова произведения квадратур с учётом того,
    // что геометрия описывается таким образом, что параметрическую зависимость имеет верхняя и нижняя границы,
    // а правая и левая заданы константами.
    void set_quadrature(const quadrature_base<Type>& quadrature_xi, const quadrature_base<Type>& quadrature_eta) override
    {
        Type jacobian_xi = (              boundary(side_2d::RIGHT, 0) -               boundary(side_2d::LEFT, 0)) /
                           (quadrature_xi.boundary(side_1d::RIGHT   ) - quadrature_xi.boundary(side_1d::LEFT   ));
        std::vector<Type> jacobian_eta(quadrature_xi.nodes_count()),
                          xi(quadrature_xi.nodes_count()),
                          eta(quadrature_eta.nodes_count());

        weights.resize(quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        for(size_t i = 0; i < quadrature_xi.nodes_count(); ++i) {
            xi[i] = boundary(side_2d::LEFT, 0) + (quadrature_xi.node(i) - quadrature_xi.boundary(side_1d::LEFT)) * jacobian_xi;
            jacobian_eta[i] = (               boundary(side_2d::UP, xi[i]) -                boundary(side_2d::DOWN, xi[i])) /
                              (quadrature_eta.boundary(side_1d::RIGHT    ) - quadrature_eta.boundary(side_1d::LEFT       ));
            for(size_t j = 0; j < quadrature_eta.nodes_count(); ++j) {
                eta[j] = boundary(side_2d::DOWN, xi[i]) + (quadrature_eta.node(j)-quadrature_eta.boundary(side_1d::LEFT)) * jacobian_eta[i];
                weights[i*quadrature_eta.nodes_count() + j] = quadrature_xi.weight(i) * jacobian_xi * quadrature_eta.weight(j) * jacobian_eta[i];
            }
        }

        N_in_quad_nodes.resize(element_2d<Type, Element_Type>::nodes_count(), quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        Nxi_in_quad_nodes.resize(element_2d<Type, Element_Type>::nodes_count(), quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        Neta_in_quad_nodes.resize(element_2d<Type, Element_Type>::nodes_count(), quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        for(size_t i = 0; i < nodes_count(); ++i)
            for(size_t j = 0; j < quadrature_xi.nodes_count(); ++j)
                for(size_t k = 0; k < quadrature_eta.nodes_count(); ++k) {
                    N_in_quad_nodes   (i, j*quadrature_eta.nodes_count() + k) = N   (i, xi[j], eta[k]);
                    Nxi_in_quad_nodes (i, j*quadrature_eta.nodes_count() + k) = Nxi (i, xi[j], eta[k]);
                    Neta_in_quad_nodes(i, j*quadrature_eta.nodes_count() + k) = Neta(i, xi[j], eta[k]);
                }
    }
};

}

#endif
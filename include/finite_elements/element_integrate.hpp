#ifndef FINITE_ELEMENT_INTEGRATE_HPP
#define FINITE_ELEMENT_INTEGRATE_HPP

#include "element_1d.hpp"
#include "element_2d.hpp"

namespace finite_element {

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

    explicit element_1d_integrate(const quadrature_base<Type> &quadrature) { set_quadrature(quadrature); }

    void set_quadrature(const quadrature_base<Type> &quadrature) override {
        std::vector<Type> xi(quadrature.nodes_count());
        weights.resize(quadrature.nodes_count());
        Type jacobian = (           boundary(side_1d::RIGHT) -            boundary(side_1d::LEFT)) /
                        (quadrature.boundary(side_1d::RIGHT) - quadrature.boundary(side_1d::LEFT));
        for(size_t q = 0; q < quadrature.nodes_count(); ++q) {
            xi[q] = boundary(side_1d::LEFT) + (quadrature.node(q) - quadrature.boundary(side_1d::LEFT)) * jacobian;
            std::cout << quadrature.weight(q) << " " << jacobian << std::endl;
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

    explicit element_2d_integrate(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) {
        set_quadrature(quadrature_xi, quadrature_eta);
    }

    // Узлы в которых происходит расчёт получаются путём декартова произведения квадратур с учётом того,
    // что геометрия описывается таким образом, что параметрическую зависимость имеет верхняя и нижняя границы,
    // а правая и левая заданы константами.
    void set_quadrature(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) override
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

/*
// Частичная специализация класса описанного выше, на случай кубических серендиповых элементов. Необходима для экспериментальных целей.
template<class Type>
class element_2d_integrate<Type, qubic_serendip> : public element_2d_integrate_base<Type>,
                                                   public element_2d<Type, qubic_serendip> {
    std::vector<Type> weights, N_in_quad_nodes, Nxi_in_quad_nodes, Neta_in_quad_nodes,
                      Nxi2InQuad, NxietaInQuad, Neta2InQuad, Nxi3InQuad, Nxi2etaInQuad, Nxieta2InQuad, Neta3InQuad;

public:
    element_2d_integrate(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) { set(quadrature_xi, quadrature_eta); }

    // Узлы в которых происходит расчёт получаются путём декартова произведения квадратур с учётом того,
    // что геометрия описывается таким образом, что параметрическую зависимость имеет верхняя и нижняя границы,
    // а правая и левая заданы константами.
    void set(const quadrature_base<Type> &quadrature_xi, const quadrature_base<Type> &quadrature_eta) override
    {
        Type jacobian_xi = (this-> boundary(side_2d::RIGHT, 0) - this ->boundary(side_2d::LEFT, 0)) /
                          (quadrature_xi.boundary(side_1d::RIGHT)    - quadrature_xi.boundary(side_1d::LEFT));
        std::vector<Type> jacobian_eta(quadrature_xi. nodes_count()),
                          xi(quadrature_xi. nodes_count()),
                          eta(quadrature_eta.nodes_count());

        weights.resize(quadrature_xi.nodes_count() * quadrature_eta.nodes_count());
        for(size_t i = 0; i < quadrature_xi.nodes_count(); ++i) {
            xi[i] = this->boundary(side_2d::LEFT, 0) + (quadrature_xi.node(i) - quadrature_xi.boundary(side_1d::LEFT)) * jacobian_xi;
            jacobian_eta[i] = (this  ->boundary(side_2d::UP, xi[i]) - this  ->boundary(side_2d::DOWN, xi[i])) /
                             (quadrature_eta.boundary(side_1d::RIGHT)     - quadrature_eta.boundary(side_1d::LEFT));
            for(size_t j = 0; j < quadrature_eta.nodes_count(); ++j) {
                eta[j] = this->boundary(side_2d::DOWN, xi[i]) + (quadrature_eta.node(j)-quadrature_eta.boundary(side_1d::LEFT)) * jacobian_eta[i];
                weights[i*quadrature_eta.nodes_count() + j] = quadrature_xi.weight(i) * jacobian_xi * quadrature_eta.weight(j) * jacobian_eta[i];
            }
        }

        N_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());
        Neta_in_quad_nodes.resize(int_base::qnodes_count() * this->nodes_count());

        Nxi2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        NxietaInQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Neta2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi3InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxi2etaInQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Nxieta2InQuad.resize(int_base::qnodes_count() * this->nodes_count());
        Neta3InQuad.resize(int_base::qnodes_count() * this->nodes_count());

        for(size_t i = 0; i < this->nodes_count(); ++i)
            for(size_t j = 0; j < quadrature_xi.nodes_count(); ++j)
                for(size_t k = 0; k < quadrature_eta.nodes_count(); ++k) {
                    N_in_quad_nodes   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->N   (i, xi[j], eta[k]);
                    Nxi_in_quad_nodes [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi (i, xi[j], eta[k]);
                    Neta_in_quad_nodes[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta(i, xi[j], eta[k]);

                    Nxi2InQuad   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi2   (i, xi[j], eta[k]);
                    NxietaInQuad [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxieta (i, xi[j], eta[k]);
                    Neta2InQuad  [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta2  (i, xi[j], eta[k]);
                    Nxi3InQuad   [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi3   (i, xi[j], eta[k]);
                    Nxi2etaInQuad[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxi2eta(i, xi[j], eta[k]);
                    Nxieta2InQuad[i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Nxieta2(i, xi[j], eta[k]);
                    Neta3InQuad  [i*int_base::qnodes_count() + j*quadrature_eta.nodes_count() + k] = this->Neta3  (i, xi[j], eta[k]);
                }
    }

    size_t int_base::qnodes_count() const override { return weights.size(); }

    Type weight(const size_t q) const override { return weights[q]; }

    Type qN   (const size_t i, const size_t q) const override { return N_in_quad_nodes   [i*int_base::qnodes_count() + q]; }
    Type qNxi (const size_t i, const size_t q) const override { return Nxi_in_quad_nodes [i*int_base::qnodes_count() + q]; }
    Type qNeta(const size_t i, const size_t q) const override { return Neta_in_quad_nodes[i*int_base::qnodes_count() + q]; }

    Type qNxi2   (const size_t i, const size_t q) const { return Nxi2InQuad   [i*int_base::qnodes_count() + q]; }
    Type qNxieta (const size_t i, const size_t q) const { return NxietaInQuad [i*int_base::qnodes_count() + q]; }
    Type qNeta2  (const size_t i, const size_t q) const { return Neta2InQuad  [i*int_base::qnodes_count() + q]; }
    Type qNxi3   (const size_t i, const size_t q) const { return Nxi3InQuad   [i*int_base::qnodes_count() + q]; }
    Type qNxi2eta(const size_t i, const size_t q) const { return Nxi2etaInQuad[i*int_base::qnodes_count() + q]; }
    Type qNxieta2(const size_t i, const size_t q) const { return Nxieta2InQuad[i*int_base::qnodes_count() + q]; }
    Type qNeta3  (const size_t i, const size_t q) const { return Neta3InQuad  [i*int_base::qnodes_count() + q]; }

    virtual ~element_2d_integrate() = default;
};
*/

}

#endif
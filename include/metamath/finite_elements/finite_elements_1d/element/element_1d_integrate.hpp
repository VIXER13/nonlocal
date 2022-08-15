#ifndef FINITE_ELEMENT_1D_INTEGRATE_HPP
#define FINITE_ELEMENT_1D_INTEGRATE_HPP

#include "element_1d.hpp"
#include "element_1d_integrate_base.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Element_Type, auto... Args>
class element_1d_integrate : public element_1d_integrate_base<T>,
                             public element_1d<T, Element_Type, Args...> {
    using element_1d_t = element_1d<T, Element_Type, Args...>;
    using element_integrate_1d_t = element_1d_integrate_base<T>;
    using element_integrate_1d_t::_quadrature;
    using element_integrate_1d_t::_weights;
    using element_integrate_1d_t::_qN;
    using element_integrate_1d_t::_qNxi;

public:
    using element_integrate_1d_t::qnodes_count;
    using element_integrate_1d_t::nodes_count;
    using element_1d_t::boundary;
    using element_1d_t::N;
    using element_1d_t::Nxi;

    explicit element_1d_integrate(const quadrature_1d_base<T>& quadrature) { set_quadrature(quadrature); }
    ~element_1d_integrate() override = default;

    void set_quadrature(const quadrature_1d_base<T>& quadrature) override {
        _quadrature = quadrature.clone();
        std::vector<T> xi(quadrature.nodes_count());
        _weights.resize(quadrature.nodes_count());
        T jacobian = (           boundary(side_1d::RIGHT) -            boundary(side_1d::LEFT)) /
                     (quadrature.boundary(side_1d::RIGHT) - quadrature.boundary(side_1d::LEFT));
        for(size_t q = 0; q < quadrature.nodes_count(); ++q) {
            xi[q] = boundary(side_1d::LEFT) + (quadrature.node(q) - quadrature.boundary(side_1d::LEFT)) * jacobian;
            _weights[q] = quadrature.weight(q) * jacobian;
        }

        _qN  .resize(element_1d_t::nodes_count() * qnodes_count());
        _qNxi.resize(element_1d_t::nodes_count() * qnodes_count());
        for(size_t i = 0; i < nodes_count(); ++i)
            for(size_t q = 0; q < qnodes_count(); ++q) {
                _qN  [i*qnodes_count() + q] = N  (i, xi[q]);
                _qNxi[i*qnodes_count() + q] = Nxi(i, xi[q]);
            }
    }
};

}

#endif
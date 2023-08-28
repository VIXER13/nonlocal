#ifndef FINITE_ELEMENT_1D_INTEGRATE_HPP
#define FINITE_ELEMENT_1D_INTEGRATE_HPP

#include "element_1d.hpp"
#include "element_1d_integrate_base.hpp"

namespace metamath::finite_element {

template<class T, size_t Derivative_Order, template<class, auto...> class Element_Type, auto... Args>
class element_1d_integrate : public element_1d_integrate_base<T>,
                             public element_1d<T, Derivative_Order, Element_Type, Args...> {
    using element_1d_t = element_1d<T, Derivative_Order, Element_Type, Args...>;
    using element_integrate_1d_t = element_1d_integrate_base<T>;
    using element_integrate_1d_t::_nearest_qnode;
    using element_integrate_1d_t::_quadrature;
    using element_integrate_1d_t::_weights;
    using element_integrate_1d_t::_qN;
    using element_integrate_1d_t::qnodes_by_nodes;

public:
    using element_integrate_1d_t::max_derivative_order;
    using element_integrate_1d_t::qnodes_count;
    using element_integrate_1d_t::nodes_count;
    using element_integrate_1d_t::qnodes;
    using element_1d_t::boundary;
    using element_1d_t::nodes;
    using element_1d_t::node;
    using element_1d_t::N;

    explicit element_1d_integrate(const quadrature_1d_base<T>& quadrature)
        : element_integrate_1d_t{Derivative_Order} {
        set_quadrature(quadrature);
    }

    ~element_1d_integrate() override = default;

    void set_quadrature(const quadrature_1d_base<T>& quadrature) override {
        using enum side_1d;
        _quadrature = quadrature.clone();
        _weights.resize(quadrature.nodes_count());
        _qN.resize((max_derivative_order() + 1) * element_1d_t::nodes_count() * qnodes_count());
        std::vector<T> qnodes_coord(quadrature.nodes_count());
        const T jacobian = (           boundary(RIGHT) -            boundary(LEFT)) /
                           (quadrature.boundary(RIGHT) - quadrature.boundary(LEFT));

        for(const size_t q : qnodes()) {
            _weights[q] = quadrature.weight(q) * jacobian;
            qnodes_coord[q] = boundary(LEFT) + (quadrature.node(q) - quadrature.boundary(LEFT)) * jacobian;
            for(const size_t d : std::ranges::iota_view{0u, max_derivative_order() + 1})
                for(const size_t i : nodes())
                    _qN[d*qnodes_by_nodes() + i*qnodes_count() + q] = N(d, i, qnodes_coord[q]);
        }

        _nearest_qnode.resize(element_1d_t::nodes_count());
        for(const size_t i : nodes()) {
            T length = std::numeric_limits<T>::max();
            for(const size_t q : qnodes())
                if (const T curr_length = std::abs(qnodes_coord[q] - node(i)); length > curr_length) {
                    length = curr_length;
                    _nearest_qnode[i] = q;
                }
        }
    }
};

}

#endif
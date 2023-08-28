#ifndef FINITE_ELEMENT_1D_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_1D_INTEGRATE_BASE_HPP

#include "finite_element_integrate_base.hpp"
#include "element_1d_base.hpp"
#include "quadrature_1d_base.hpp"

namespace metamath::finite_element {

template<std::floating_point T>
class element_integrate_base_1d {
protected:
    std::vector<size_t> _nearest_qnode;
    std::vector<T> _weights;
    std::vector<T> _qN;
    const size_t _max_derivative_order;

    explicit element_integrate_base_1d(const size_t max_derivative_order) noexcept
        : _max_derivative_order{max_derivative_order} {}

    size_t qnodes_by_nodes() const noexcept {
        return qnodes_count() * nodes_count();
    }

public:
    virtual ~element_integrate_base_1d() noexcept = default;

    size_t max_derivative_order() const noexcept {
        return _max_derivative_order;
    }

    size_t qnodes_count() const noexcept {
        return _weights.size();
    }

    size_t  nodes_count() const noexcept {
        return _qN.size() / ((max_derivative_order() + 1) * qnodes_count());
    }

    std::ranges::iota_view<size_t, size_t> qnodes() const noexcept {
        return {0u, qnodes_count()};
    }

    std::ranges::iota_view<size_t, size_t>  nodes() const noexcept {
        return {0u,  nodes_count()};
    }

    size_t nearest_qnode(const size_t i) const {
        return _nearest_qnode[i];
    }

    T weight(const size_t q) const {
        return _weights[q];
    }

    T qN(const size_t d, const size_t i, const size_t q) const {
        return _qN[d*qnodes_by_nodes() + i*qnodes_count() + q];
    }
};

template<std::floating_point T>
class element_1d_integrate_base : public element_integrate_base_1d<T>,
                                  private virtual element_1d_base<T> {
protected:
    std::unique_ptr<quadrature_1d_base<T>> _quadrature = nullptr;

    explicit element_1d_integrate_base(const size_t max_derivative_order) noexcept
        : element_integrate_base_1d<T>{max_derivative_order} {}

public:
    using element_integrate_base_1d<T>::max_derivative_order;
    using element_integrate_base_1d<T>::qnodes_count;
    using element_integrate_base_1d<T>::nodes_count;
    using element_integrate_base_1d<T>::nodes;
    using element_integrate_base_1d<T>::qnodes;
    using element_1d_base<T>::boundary;
    using element_1d_base<T>::node;
    using element_1d_base<T>::N;

    ~element_1d_integrate_base() override = default;

    virtual void set_quadrature(const quadrature_1d_base<T>& quadrature) = 0;

    const quadrature_1d_base<T>& quadrature() const {
        return *_quadrature;
    }
};

}

#endif
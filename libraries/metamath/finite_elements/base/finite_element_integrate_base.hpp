#pragma once

#include <vector>

namespace metamath::finite_element {

template<class T>
class element_integrate_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    std::vector<T> _weights, _qN;
    std::vector<size_t> _nearest_qnode;

    explicit element_integrate_base() noexcept = default;

public:
    virtual ~element_integrate_base() noexcept = default;

    size_t qnodes_count() const noexcept { return _weights.size(); }
    size_t  nodes_count() const noexcept { return _qN.size() / qnodes_count(); }
    std::ranges::iota_view<size_t, size_t> qnodes() const noexcept { return {0u, qnodes_count()}; }
    std::ranges::iota_view<size_t, size_t>  nodes() const noexcept { return {0u,  nodes_count()};}

    size_t nearest_qnode(const size_t i) const { return _nearest_qnode[i]; }
    T weight(const size_t q) const { return _weights[q]; }
    T qN(const size_t i, const size_t q) const { return _qN[i*qnodes_count() + q]; }
};

}
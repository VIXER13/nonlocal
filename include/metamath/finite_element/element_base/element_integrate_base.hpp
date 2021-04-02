#ifndef FINITE_ELEMENT_INTEGRATE_BASE_HPP
#define FINITE_ELEMENT_INTEGRATE_BASE_HPP

#include <vector>

namespace metamath::finite_element {

template<class T>
class element_integrate_base {
    static_assert(std::is_floating_point_v<T>, "The T must be floating point.");

protected:
    std::vector<T> _weights, _qN;
    explicit element_integrate_base() noexcept = default;

public:
    virtual ~element_integrate_base() noexcept = default;

    size_t qnodes_count() const noexcept { return _weights.size(); }
    size_t  nodes_count() const noexcept { return _qN.size() / qnodes_count(); }

    T weight(const size_t q) const noexcept { return _weights[q]; }

    T qN(const size_t i, const size_t q) const noexcept { return _qN[i*qnodes_count() + q]; }
};

}

#endif
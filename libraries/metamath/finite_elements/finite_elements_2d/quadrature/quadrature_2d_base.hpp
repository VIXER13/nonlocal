#pragma once

#include <metamath/finite_elements/base/quadrature_base.hpp>
#include <metamath/finite_elements/finite_elements_2d/geometry/side_2d.hpp>

#include <memory>

namespace metamath::finite_element {

template<std::floating_point T>
class quadrature_2d_base : public quadrature_base<T> {
public:
    ~quadrature_2d_base() override = default;
    virtual std::unique_ptr<quadrature_2d_base<T>> copy() const = 0;
    virtual const std::array<T, 2>& node(const size_t i) const = 0;
    virtual T boundary(const side_2d bound, const T x) const = 0;
};

}
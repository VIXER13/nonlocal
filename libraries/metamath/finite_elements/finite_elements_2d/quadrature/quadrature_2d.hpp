#pragma once

#include "quadrature_2d_base.hpp"

namespace metamath::finite_element {

template<std::floating_point T, template<class, auto...> class Quadrature_Type, auto... Args>
class quadrature_2d : public quadrature_2d_base<T>,
                      public Quadrature_Type<T, Args...> {
    using quadrature_t = Quadrature_Type<T, Args...>;
    static_assert(quadrature_t::nodes.size() == quadrature_t::weights.size(),
                  "Number of nodes must match number of weights");

public:
    ~quadrature_2d() override = default;

    std::unique_ptr<quadrature_2d_base<T>> copy() const override {
        return std::make_unique<quadrature_2d<T, Quadrature_Type, Args...>>(*this);
    }

    size_t nodes_count() const override { return quadrature_t::nodes.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return quadrature_t::nodes[i]; }
    T weight(const size_t i) const override { return quadrature_t::weights[i]; }
    T boundary(const side_2d bound, const T x) const override { return quadrature_t::boundary(bound, x); }
};

}
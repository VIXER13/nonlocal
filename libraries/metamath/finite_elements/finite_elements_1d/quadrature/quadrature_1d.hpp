#pragma once

#include "quadrature_1d_base.hpp"

namespace metamath::finite_element {

template<class T, template<class, auto...> class Quadrature_Type, auto... Args>
class quadrature_1d : public quadrature_1d_base<T>,
                      public Quadrature_Type<T, Args...> {
    using quadrature_t = Quadrature_Type<T, Args...>;
    static_assert(quadrature_t::nodes.size() == quadrature_t::weights.size(),
                  "The number of nodes and weights does not match.");

public:
    ~quadrature_1d() override = default;

    std::unique_ptr<quadrature_1d_base<T>> clone() const override { return std::make_unique<quadrature_1d<T, Quadrature_Type, Args...>>(); }

    size_t nodes_count() const override { return quadrature_t::nodes.size(); }

    const T node(const size_t i) const override { return quadrature_t::nodes[i]; }
    T weight(const size_t i) const override { return quadrature_t::weights[i]; }

    T boundary(const side_1d bound) const override { return quadrature_t::boundary(bound); }
};

}
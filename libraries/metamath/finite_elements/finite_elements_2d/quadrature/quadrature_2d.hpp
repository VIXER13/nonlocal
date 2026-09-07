#pragma once

#include "quadrature_2d_base.hpp"
#include "primitives/cartesian_production.hpp"

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

template<std::floating_point T>
class quadrature_2d<T, cartesian_production> : public quadrature_2d_base<T>
                                             , public cartesian_production<T> {
    using quadrature_t = cartesian_production<T>;

public:
    explicit quadrature_2d(const quadrature_1d_base<T>& quadrature)
        : quadrature_2d<T, cartesian_production>{quadrature, quadrature} {}
    explicit quadrature_2d(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y)
        : cartesian_production<T>{quadrature_x, quadrature_y} {
        if (quadrature_t::nodes.size() != quadrature_t::weights.size())
            throw std::runtime_error{"Number of nodes shall match number of weights"};
    }
    ~quadrature_2d() override = default;

    std::unique_ptr<quadrature_2d_base<T>> copy() const override {
        return std::make_unique<quadrature_2d<T, cartesian_production>>(*this);
    }

    size_t nodes_count() const override { return quadrature_t::nodes.size(); }

    const std::array<T, 2>& node(const size_t i) const override { return quadrature_t::nodes[i]; }
    T weight(const size_t i) const override { return quadrature_t::weights[i]; }
    T boundary(const side_2d bound, const T x) const override { return quadrature_t::boundary(bound, x); }
};

}
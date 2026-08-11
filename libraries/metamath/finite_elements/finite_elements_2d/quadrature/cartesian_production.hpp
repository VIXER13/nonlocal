#pragma once

#include "quadrature_2d.hpp"

#include <metamath/finite_elements/finite_elements_1d/quadrature/quadrature_1d.hpp>
#include <metamath/finite_elements/finite_elements_2d/geometry/geometry_2d.hpp>

#include <vector>
#include <ranges>

namespace metamath::finite_element {

template<std::floating_point T>
class cartesian_production : public geometry_2d<T, rectangle_element_geometry> {
protected:
    std::vector<T> weights;
    std::vector<std::array<T, 2>> nodes;

    explicit cartesian_production(const quadrature_1d_base<T>& quadrature_x, const quadrature_1d_base<T>& quadrature_y) 
        : weights(quadrature_x.nodes_count() * quadrature_y.nodes_count())
        , nodes(quadrature_x.nodes_count() * quadrature_y.nodes_count()) {
        for(const size_t i : std::ranges::iota_view{0zu, quadrature_x.nodes_count()})
            for(const size_t j : std::ranges::iota_view{0zu, quadrature_y.nodes_count()}) {
                weights[i * quadrature_y.nodes_count() + j] = quadrature_x.weight(i) * quadrature_y.weight(j);
                nodes[i * quadrature_y.nodes_count() + j] = {quadrature_x.node(i), quadrature_y.node(j)};
            }
    }
    ~cartesian_production() noexcept override = default;
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
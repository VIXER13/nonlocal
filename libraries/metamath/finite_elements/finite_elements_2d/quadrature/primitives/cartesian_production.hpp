#pragma once

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

}
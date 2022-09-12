#ifndef FINITE_ELEMENTS_BASIS_LAGRANGIAN_ELEMENTS_2D_HPP
#define FINITE_ELEMENTS_BASIS_LAGRANGIAN_ELEMENTS_2D_HPP

#include "geometry_2d.hpp"
#include "array_cartesian_product.hpp"
#include "uniform_partition.hpp"

namespace metamath::finite_element {

template<class T, size_t N, size_t M>
class lagrangian_element_2d : public geometry_2d<T, rectangle_element_geometry> {
    static inline constexpr std::array<T, N + 1> nodes_x = utils::uniform_partition<N + 1>(std::array{T{-1}, T{1}});
    static inline constexpr std::array<T, M + 1> nodes_y = utils::uniform_partition<M + 1>(std::array{T{-1}, T{1}});

protected:
    using geometry_2d<T, rectangle_element_geometry>::x;
    using geometry_2d<T, rectangle_element_geometry>::y;

    static inline constexpr std::array<std::array<T, 2>, (N + 1) * (M + 1)>
        nodes = utils::array_cartesian_product(nodes_x, nodes_y);
    static inline constexpr auto
        basis = symbolic::basis_production(
            symbolic::generate_lagrangian_basis<x>(nodes_x),
            symbolic::generate_lagrangian_basis<y>(nodes_y)
        );

    constexpr explicit lagrangian_element_2d() noexcept = default;
    ~lagrangian_element_2d() override = default;
};

}

#endif
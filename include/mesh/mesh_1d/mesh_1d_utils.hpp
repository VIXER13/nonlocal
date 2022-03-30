#ifndef MESH_1D_UTILS_HPP
#define MESH_1D_UTILS_HPP

#include "mesh_1d.hpp"
#include <ranges>

namespace nonlocal::utils {

template<class T, class Vector>
T integrate_solution(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"nodes_count() != x.size()"};
    T integral = 0;
    const auto& el = mesh.element();
    for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count()))
        for(const size_t i : std::views::iota(size_t{0}, el.nodes_count()))
            for(const size_t q : std::views::iota(size_t{0}, el.qnodes_count()))
                integral += el.weight(q) * el.qN(i, q) * x[mesh.node_number(e, i)];
    return integral * mesh.jacobian();
}

template<class T, class Vector>
std::vector<T> gradient(const mesh::mesh_1d<T>& mesh, const Vector& x) {
    if (mesh.nodes_count() != x.size())
        throw std::logic_error{"nodes_count() != x.size()"};
    std::vector dx(x.size(), T{0});
    for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count()))
        for(const size_t i : std::views::iota(size_t{0}, mesh.element().nodes_count() - (e != mesh.elements_count()-1))) {
            const size_t node = mesh.node_number(e, i);
            const T coord = mesh.node_coord(node);
            for(const size_t j : std::views::iota(size_t{0}, mesh.element().nodes_count()))
                dx[node] += mesh.element().Nxi(j, coord) * x[mesh.node_number(e, j)] / mesh.jacobian();
        }
    return dx;
}

}

#endif
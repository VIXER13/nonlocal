#ifndef NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_FIRST_KIND_2D_HPP

#include "boundary_condition_2d.hpp"
#include "mesh_2d.hpp"
#include <eigen3/Eigen/Dense>

namespace nonlocal
{

template<class T, class I, class Callback>
void boundary_nodes_run(const mesh::mesh_2d<T, I>& mesh, const Callback& callback) {
    for(const std::string& b : mesh.boundary_names())
        for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(b)))
            for(const size_t i : std::views::iota(size_t{0}, mesh.element_1d(b, e)->nodes_count()))
                callback(b, e, i);
}

template<class B, class T, class I, class Matrix_Index, size_t DoF>
void boundary_condition_first_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                      const mesh::mesh_proxy<T, I>& mesh_proxy,
                                      const std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>>& boundary_condition,
                                      const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(K_bound.cols());
    boundary_nodes_run(mesh_proxy.mesh(),
        [&mesh = mesh_proxy.mesh(), &boundary_condition, &x](const std::string& b, const size_t el, const size_t i) {
            for(const size_t comp : std::views::iota(size_t{0}, DoF))
                if (boundary_condition.at(b)[comp].type == B(boundary_condition_t::FIRST_KIND))
                    if (const I node = DoF * mesh.node_number(b, el, i) + comp; x[node] == 0)
                        x[node] = boundary_condition.at(b)[comp].func(mesh.node(node));
        });

    f -= K_bound * x;

    using namespace metamath::functions;
    const std::array<size_t, 2> first_last = {DoF * mesh_proxy.first_node(), DoF * mesh_proxy.last_node()};
    boundary_nodes_run(mesh_proxy.mesh(),
        [&mesh = mesh_proxy.mesh(), &boundary_condition, &x, &f, &first_last](const std::string& b, const size_t el, const size_t i) {
            for(const size_t comp : std::views::iota(size_t{0}, DoF))
                if (boundary_condition.at(b)[comp].type == B(boundary_condition_t::FIRST_KIND))
                    if (const I node = DoF * mesh.node_number(b, el, i) + comp;
                        node >= first_last.front() && node < first_last.back())
                        f[node - first_last.front()] = x[node];
        });
}
    
}

#endif
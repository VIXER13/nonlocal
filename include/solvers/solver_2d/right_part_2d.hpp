#ifndef NONLOCAL_RIGHT_PART_2D_HPP
#define NONLOCAL_RIGHT_PART_2D_HPP

#include "../solvers_constants.hpp"
#include "boundary_condition_2d.hpp"
#include "mesh_2d.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace nonlocal {

template<class T, class I, class Callback>
void boundary_nodes_run(const mesh::mesh_2d<T, I>& mesh, const Callback& callback) {
    for(const std::string& b : mesh.boundary_names())
        for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(b)))
            for(const size_t i : std::views::iota(size_t{0}, mesh.element_1d(b, e)->nodes_count()))
                callback(b, e, i);
}

template<class T, class I, class B, size_t DoF>
std::vector<bool> inner_nodes(const mesh::mesh_2d<T, I>& mesh,
                              const std::unordered_map<std::string, std::array<B, DoF>>& boundary_condition) {
    std::vector<bool> is_inner(DoF * mesh.nodes_count(), true);
    boundary_nodes_run(mesh,
        [&mesh, &boundary_condition, &is_inner](const std::string& b, const size_t el, const size_t i) {
            for(const size_t comp : std::views::iota(size_t{0}, DoF))
                if(boundary_condition.at(b)[comp] == B(boundary_condition_t::FIRST_KIND))
                    is_inner[DoF * mesh.node_number(b, el, i) + comp] = false;
        });
    return is_inner;
}

template<class B, class T, class I, size_t DoF>
void boundary_condition_first_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                      const mesh::mesh_proxy<T, I>& mesh_proxy,
                                      const std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>>& boundary_condition,
                                      const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(K_bound.cols());
    boundary_nodes_run(mesh_proxy.mesh(),
        [&mesh = mesh_proxy.mesh(), &boundary_condition, &x](const std::string& b, const size_t el, const size_t i) {
            for(const size_t comp : std::views::iota(size_t{0}, DoF))
                if (boundary_condition.at(b)[comp].type == B(boundary_condition_t::FIRST_KIND))
                    if (const I node = DoF * mesh.node_number(b, el, i) + comp; x[node] == 0)
                        x[node] = boundary_condition.at(b)[comp].func(mesh.node(node));
        });

    f -= K_bound * x;

    using namespace metamath::function;
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

template<class T, class I, class Boundary_Gradient>
T integrate_boundary_gradient(const mesh::mesh_proxy<T, I>& mesh_proxy,
                              const std::string& b, const size_t e, const size_t i,
                              const Boundary_Gradient& boundary_gradient) {
    T integral = 0;
    const auto& be = mesh_proxy.mesh().element_1d(b, e);
    auto qcoord = mesh_proxy.quad_coord(b, e);
    auto J      = mesh_proxy.jacobi_matrix(b, e);
    for(size_t q = 0; q < be->qnodes_count(); ++q, ++qcoord, ++J)
        integral += be->weight(q) * be->qN(i, q) * boundary_gradient(*qcoord) * mesh_proxy.jacobian(*J);
    return integral;
}

template<class B, class T, class I, size_t DoF>
void boundary_condition_second_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const mesh::mesh_proxy<T, I>& mesh_proxy,
                                       const std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>>& boundary_condition) {
    const auto& mesh = mesh_proxy.mesh();
    const std::array<size_t, 2>& first_last = {DoF * mesh_proxy.first_node(), DoF * mesh_proxy.last_node()};
    for(const std::string& b : mesh.boundary_names())
        for(const size_t comp : std::views::iota(size_t{0}, DoF))
            if(boundary_condition.at(b)[comp].type == B(boundary_condition_t::SECOND_KIND))
                for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(b)))
                    for(const size_t i : std::views::iota(size_t{0}, mesh.element_1d(b, e)->nodes_count()))
                        if(const I node = DoF * mesh.node_number(b, e, i) + comp;
                           node >= first_last.front() && node < first_last.back())
                            f[node - first_last.front()] += integrate_boundary_gradient(mesh_proxy, b, e, i, boundary_condition.at(b)[comp].func);
}

template<size_t DoF, class T, class I, class Right_Part>
void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const mesh::mesh_proxy<T, I>& mesh_proxy, const Right_Part& right_part) {
    const auto integrate_function = [&mesh_proxy, &func = right_part](const size_t e, const size_t i) {
        std::conditional_t<DoF == 1, T, std::array<T, DoF>> integral = {};
        const auto& el = mesh_proxy.mesh().element_2d(e);
        auto J = mesh_proxy.jacobi_matrix(e);
        auto qcoord = mesh_proxy.quad_coord(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J) {
            using namespace metamath::function;
            integral += el->weight(q) * el->qN(i, q) * mesh_proxy.jacobian(*J) * func(*qcoord);
        }
        return integral;
    };

#pragma omp parallel for default(none) shared(f, mesh_proxy, integrate_function)
    for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node)
        for(const I e : mesh_proxy.nodes_elements_map(node)) {
            const size_t i = mesh_proxy.global_to_local_numbering(e).find(node)->second;
            const std::array<T, DoF> integral = {integrate_function(e, i)};
            for(const size_t comp : std::ranges::iota_view(size_t{0}, DoF))
                f[DoF * (node - mesh_proxy.first_node()) + comp] += integral[comp];
        }
}

}

#endif
#ifndef CONVECTION_CONDITION_2D_HPP
#define CONVECTION_CONDITION_2D_HPP

#include "mesh_2d.hpp"
#include "right_part_2d.hpp"
#include <eigen3/Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I, class Callback>
void convection_boundary_run(const mesh::mesh_proxy<T, I>& mesh_proxy,
                             const std::unordered_map<std::string, std::array<boundary_condition_t, 1>>& bounds_types,
                             const std::unordered_map<std::string, T>& heat_transfer,
                             const Callback& callback) {
    const auto& mesh = mesh_proxy.mesh();
    for(const auto& [bound_name, condition_type] : bounds_types)
        if (condition_type.front() == boundary_condition_t::CONVECTION) {
            const T alpha = heat_transfer.at(bound_name);
            for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(bound_name)))
                callback(alpha, bound_name, e);
        }
}

template<class T, class I, class Matrix_Index>
void convection_condition_matrix_part_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                                         const mesh::mesh_proxy<T, I>& mesh_proxy,
                                         const std::unordered_map<std::string, std::array<boundary_condition_t, 1>>& boundary_condition,
                                         const std::unordered_map<std::string, T>& alpha) {
    const auto integrate = [&mesh_proxy](const std::string& b, const size_t e, const size_t i, const size_t j) {
        T integral = 0;
        auto J = mesh_proxy.jacobi_matrix(b, e);
        const auto& be = mesh_proxy.mesh().element_1d(b, e);
        for(size_t q = 0; q < be->qnodes_count(); ++q, ++J)
            integral += be->weight(q) * be->qN(i, q) * be->qN(j, q) * mesh::jacobian(*J);
        return integral;
    };

    convection_boundary_run(mesh_proxy, boundary_condition, alpha,
        [&K, &mesh_proxy, &integrate](const T alpha, const std::string& b, const size_t e) {
            const auto& mesh = mesh_proxy.mesh();
            const auto& be = mesh.element_1d(b, e);
            for(const size_t i : std::views::iota(size_t{0}, be->nodes_count()))
                for(const size_t j : std::views::iota(size_t{0}, be->nodes_count())) {
                    const size_t row = mesh.node_number(b, e, i),
                                 col = mesh.node_number(b, e, j);
                    if (col >= row && row >= mesh_proxy.first_node() && row < mesh_proxy.last_node())
                        K.coeffRef(row, col) -= alpha * integrate(b, e, i, j);
                }
        }
    );
}

template<class T, class I>
void convection_condition_right_part_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                        const mesh::mesh_proxy<T, I>& mesh_proxy,
                                        const std::unordered_map<std::string, stationary_boundary_2d_t<boundary_condition_t, T, 1>>& boundary_condition,
                                        const std::unordered_map<std::string, T>& heat_transfer) {
    convection_boundary_run(mesh_proxy, boundary_type(boundary_condition), heat_transfer,
        [&f, &mesh_proxy, &boundary_condition](const T alpha, const std::string& b, const size_t e) {
            const auto& mesh = mesh_proxy.mesh();
            const auto& be = mesh.element_1d(b, e);
            for(const size_t i : std::views::iota(size_t{0}, be->nodes_count())) {
                const size_t row = mesh.node_number(b, e, i);
                if (row >= mesh_proxy.first_node() && row < mesh_proxy.last_node())
                    f[row - mesh_proxy.first_node()] += alpha * integrate_boundary_gradient(mesh_proxy, b, e, i, boundary_condition.at(b)[0].func);
            }
        }
    );
}

}

#endif
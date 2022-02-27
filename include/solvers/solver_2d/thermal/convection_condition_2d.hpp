#ifndef CONVECTION_CONDITION_2D_HPP
#define CONVECTION_CONDITION_2D_HPP

#include "mesh_2d.hpp"
#include "boundary_condition_2d.hpp"
#include <eigen3/Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
void convection_condition_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                             const mesh::mesh_proxy<T, I>& mesh_proxy,
                             const std::unordered_map<std::string, std::array<boundary_condition_t, 1>>& bounds_types,
                             const std::unordered_map<std::string, T>& alphas) {
    const auto integrate = [&mesh_proxy](const std::string& b, const size_t e, const size_t i, const size_t j) {
        T integral = 0;
        auto J = mesh_proxy.jacobi_matrix(b, e);
        const auto& be = mesh_proxy.mesh().element_1d(b, e);
        for(size_t q = 0; q < be->qnodes_count(); ++q, ++J)
            integral += be->weight(q) * be->qN(i, q) * be->qN(j, q) * mesh_proxy.jacobian(*J);
        return integral;
    };

    const auto& mesh = mesh_proxy.mesh();
    for(const auto& [bound_name, condition_type] : bounds_types)
        if (condition_type.front() == boundary_condition_t::CONVECTION) {
            const T alpha = alphas.at(bound_name);
            for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(bound_name))) {
                const auto& be = mesh.element_1d(bound_name, e);
                for(const size_t i : std::views::iota(size_t{0}, be->nodes_count()))
                    for(const size_t j : std::views::iota(size_t{0}, be->nodes_count())) {
                        const size_t row = mesh.node_number(bound_name, e, i),
                                     col = mesh.node_number(bound_name, e, j);
                        if (col >= row && row >= mesh_proxy.first_node() && row < mesh_proxy.last_node())
                            K.coeffRef(row, col) += alpha * integrate(bound_name, e, i, j);
                    }
            }
        }
}

}

#endif
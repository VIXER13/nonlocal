#pragma once

#include "thermal_boundary_conditions_2d.hpp"

#include <solvers/solver_2d/base/solvers_utils.hpp>

#include <Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
void convection_condition_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                             const mesh::mesh_2d<T, I>& mesh,
                             const thermal_boundaries_conditions_2d<T>& boundaries_conditions) {
    static constexpr auto integrate = [](const convection_2d<T>& condition, const auto& element, const size_t i, const size_t j) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * el.qN(j, q) * mesh::jacobian(element.jacobi_matrix(q));
        return condition.heat_transfer() * integral;
    };

    utils::run_by_boundaries<convection_2d>(mesh.container(), boundaries_conditions,
        [&K, &mesh, process_nodes = mesh.process_nodes()](const convection_2d<T>& condition, const size_t be, const size_t row, const size_t) {
            if (row >= process_nodes.front() && row <= process_nodes.back())
                for(const size_t j : std::ranges::iota_view{0u, mesh.container().nodes_count(be)})
                    if (const size_t col = mesh.container().node_number(be, j); col >= row)
                        K.coeffRef(row, col) += integrate(condition, mesh.container().element_1d_data(be), mesh.global_to_local(be, row), j);
        });
}

}
#pragma once

#include "thermal_boundary_conditions_2d.hpp"

#include <solvers/base/utils.hpp>

#include <Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I, class Matrix_Index>
void radiation_condition_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                            const mesh::mesh_2d<T, I>& mesh,
                            const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev,
                            const T time_step) {
    const auto integrate_matrix = [&temperature_prev](const radiation_2d<T>& condition, const auto& element, const size_t i, const size_t j) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        using namespace metamath::functions;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * power<3>(element.approximate_in_qnode(q, temperature_prev)) *
                        el.qN(i, q) * el.qN(j, q) * mesh::jacobian(element.jacobi_matrix(q));
        static constexpr T Stefan_Boltzmann_Constant_X4 = 4 * Stefan_Boltzmann_Constant<T>;
        return Stefan_Boltzmann_Constant_X4 * condition.emissivity() * integral;
    };

    const auto integrate_vector = [&temperature_prev](const radiation_2d<T>& condition, const auto& element, const size_t i) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        using namespace metamath::functions;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * power<4>(element.approximate_in_qnode(q, temperature_prev)) *
                        mesh::jacobian(element.jacobi_matrix(q));
        static constexpr T Stefan_Boltzmann_Constant_X3 = 3 * Stefan_Boltzmann_Constant<T>;
        return Stefan_Boltzmann_Constant_X3 * condition.emissivity() * integral;
    };

    utils::run_by_boundaries<radiation_2d>(mesh.container(), boundaries_conditions,
        [&K, &f, &mesh, &integrate_matrix, &integrate_vector, time_step, process_nodes = mesh.process_nodes()]
        (const radiation_2d<T>& condition, const size_t be, const size_t row, const size_t) {
            if (row >= process_nodes.front() && row <= process_nodes.back()) {
                const auto element = mesh.container().element_1d_data(be);
                const size_t i = mesh.global_to_local(be, row);
                for(const size_t j : std::ranges::iota_view{0u, mesh.container().nodes_count(be)})
                    if (const size_t col = mesh.container().node_number(be, j); col >= row)
                        K.coeffRef(row, col) += time_step * integrate_matrix(condition, element, i, j);
                f[row] += integrate_vector(condition, element, i);
            }
        });
}


template<class T, class I, class Matrix_Index>
void radiation_condition_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& residual,
                            const mesh::mesh_2d<T, I>& mesh,
                            const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev) {
    const auto integrate_matrix = [&temperature_prev](const radiation_2d<T>& condition, const auto& element, const size_t i, const size_t j) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        using namespace metamath::functions;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * power<3>(element.approximate_in_qnode(q, temperature_prev)) *
                        el.qN(i, q) * el.qN(j, q) * mesh::jacobian(element.jacobi_matrix(q));
        static constexpr T Stefan_Boltzmann_Constant_X4 = 4 * Stefan_Boltzmann_Constant<T>;
        return Stefan_Boltzmann_Constant_X4 * condition.emissivity() * integral;
    };

    const auto integrate_vector = [&temperature_prev](const radiation_2d<T>& condition, const auto& element, const size_t i) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        using namespace metamath::functions;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * power<4>(element.approximate_in_qnode(q, temperature_prev)) *
                        mesh::jacobian(element.jacobi_matrix(q));
        return Stefan_Boltzmann_Constant<T> * condition.emissivity() * integral;
    };

    utils::run_by_boundaries<radiation_2d>(mesh.container(), boundaries_conditions,
        [&K, &residual, &mesh, &integrate_matrix, &integrate_vector, process_nodes = mesh.process_nodes()]
        (const radiation_2d<T>& condition, const size_t be, const size_t row, const size_t) {
            if (row >= process_nodes.front() && row <= process_nodes.back()) {
                const auto element = mesh.container().element_1d_data(be);
                const size_t i = mesh.global_to_local(be, row);
                for(const size_t j : std::ranges::iota_view{0u, mesh.container().nodes_count(be)})
                    if (const size_t col = mesh.container().node_number(be, j); col >= row) 
                        K.coeffRef(row, col) -= integrate_matrix(condition, element, i, j);
                residual[row] -= integrate_vector(condition, element, i);
            }
        });

}

}
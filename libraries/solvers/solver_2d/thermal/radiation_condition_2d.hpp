#ifndef NONLOCAL_RADIATION_CONDITION_2D_HPP
#define NONLOCAL_RADIATION_CONDITION_2D_HPP

#include "thermal_boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include <eigen3/Eigen/Sparse>

namespace nonlocal::thermal {

template<class T, class I>
T temperature_appr(const mesh::mesh_container_2d<T, I>& mesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& temp,
                   const auto& el ,const size_t be, const size_t i) {
    T approximation = T{0};
    for(size_t q = 0; q < el.qnodes_count(); ++q){
            approximation += temp[mesh.node_number(be, i)] * el.qN(i, q);
    }
    return approximation;
}

template<class T, class I>
T radiation_rhs(const mesh::mesh_container_2d<T, I>& mesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& temp, 
                const auto& el, const T emissivity,
                const size_t be, const size_t i) {
    static constexpr T STEFAN_BOLTZMANN_CONSTANT_X3 = T{3} * STEFAN_BOLTZMANN_CONSTANT<T>;
    return STEFAN_BOLTZMANN_CONSTANT_X3 * emissivity * metamath::functions::power<4>(temperature_appr(mesh, temp, el, be, i));
}

template<class T, class I>
T radiation_mtr(const mesh::mesh_container_2d<T, I>& mesh, const Eigen::Matrix<T, Eigen::Dynamic, 1>& temp,
                const auto& el,  const T emissivity,
                const size_t be, const size_t i) {
    static constexpr T STEFAN_BOLTZMANN_CONSTANT_X4 = T{4} * STEFAN_BOLTZMANN_CONSTANT<T>;
    return STEFAN_BOLTZMANN_CONSTANT_X4 * emissivity * metamath::functions::power<3>(temperature_appr(mesh, temp, el, be, i));
}



template<class T, class I, class Matrix_Index>
void radiation_condition_2d(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                            Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                            const mesh::mesh_2d<T, I>& mesh,
                            const thermal_boundaries_conditions_2d<T>& boundaries_conditions,
                            const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature_prev,
                            const T time_step) {
    const auto integrate_matrix = [&mesh = mesh.container(), &temperature_prev](const radiation_2d<T>& condition, const size_t be, const size_t i, const size_t j) {
        T integral = T{0};
        const auto el_data = mesh.element_1d_data(be);
        const auto& [_, __, el] = el_data;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * radiation_mtr(mesh, temperature_prev, el, condition.emissivity(), be, i) * 
                        el.qN(i, q) * el.qN(j, q) * mesh::jacobian(el_data.jacobi_matrix(q));
        return integral;
    };

    const auto integrate_rhs = [&mesh = mesh.container(), &temperature_prev](const radiation_2d<T>& condition, const size_t be, const size_t i) {
        T integral = T{0};
        const auto el_data = mesh.element_1d_data(be);
        const auto& [_, __, el] = el_data;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * radiation_rhs(mesh, temperature_prev, el, condition.emissivity(), be, i) * mesh::jacobian(el_data.jacobi_matrix(q));
        return integral;
    };

    utils::run_by_boundaries<radiation_2d>(mesh.container(), boundaries_conditions,
        [&K, &f, &mesh, &integrate_matrix, &integrate_rhs, time_step, process_nodes = mesh.process_nodes()]
        (const radiation_2d<T>& condition, const size_t be, const size_t row, const size_t) {
            if (row >= process_nodes.front() && row <= process_nodes.back()){
                for(const size_t j : std::ranges::iota_view{0u, mesh.container().nodes_count(be)})
                    if (const size_t col = mesh.container().node_number(be, j); col >= row)
                        K.coeffRef(row, col) += time_step * integrate_matrix(condition, be, mesh.global_to_local(be, row), j);
                f[row] += integrate_rhs(condition, be, mesh.global_to_local(be, row));
            }
        });

}

}

#endif
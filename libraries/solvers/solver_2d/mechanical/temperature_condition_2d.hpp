#ifndef NONLOCAL_TEMPERATURE_CONDITION_2D_HPP
#define NONLOCAL_TEMPERATURE_CONDITION_2D_HPP

#include "mesh_2d.hpp"
#include "mechanical_solution_2d.hpp"
#include "../solvers_utils.hpp"
#include <eigen3/Eigen/Dense>

namespace nonlocal::mechanical {

template<class T, class I, class Influence_Function>
void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                           const mesh::mesh_proxy<T, I>& mesh_proxy,
                           const mechanical::equation_parameters<T>& parameters,
                           const T p1,
                           const Influence_Function& influence_function) {
    const T factor = parameters.thermal_expansion * parameters.E() / (T{1} - parameters.nu());
    using namespace metamath::functions;
    const std::vector<T> temperature_in_qnodes = factor * mesh::approximate_in_qnodes(mesh_proxy, parameters.delta_temperature);

    const auto integrate_temperature_loc = [&mesh_proxy, &temperature_in_qnodes](const size_t e, const size_t i) {
        std::array<T, 2> integral = {};
        auto dNd = mesh_proxy.dNdX(e, i);
        auto J = mesh_proxy.jacobi_matrix(e);
        const auto& el = mesh_proxy.mesh().element_2d(e);
        for(size_t q = 0, qshift = mesh_proxy.quad_shift(e); q < el->qnodes_count(); ++q, ++dNd, ++qshift, ++J)
            for(const size_t comp : std::ranges::iota_view{0, 2})
                integral[comp] += el->weight(q) * (*dNd)[comp] * temperature_in_qnodes[qshift];
        return integral;
    };

#pragma omp parallel for default(none) shared(f, p1, mesh_proxy, integrate_temperature_loc)
    for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node) {
        std::array<T, 2> integral = {};
        for(const I e : mesh_proxy.nodes_elements_map(node))
            integral += integrate_temperature_loc(e, mesh_proxy.global_to_local_numbering(e, node));
        for(const size_t comp : std::ranges::iota_view{0, 2})
            f[2 * node + comp] += integral[comp];
    }

    if (p1 < MAX_NONLOCAL_WEIGHT<T>) {
        const auto integrate_temperature_nonloc = [&mesh_proxy, &temperature_in_qnodes, &influence_function](const size_t eL, const size_t eNL, const size_t iL) {
            std::array<T, 2> integral = {};
            auto dNdL = mesh_proxy.dNdX(eL, iL);
            auto qcoordL = mesh_proxy.quad_coord(eL);
            const auto& elL = mesh_proxy.mesh().element_2d(eL);
            const auto& elNL = mesh_proxy.mesh().element_2d(eNL);
            const auto JNL_begin = mesh_proxy.jacobi_matrix(eNL);
            const auto qcoordNL_begin = mesh_proxy.quad_coord(eNL);
            const size_t qshiftNL_begin = mesh_proxy.quad_shift(eNL);
            for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
                auto JNL = JNL_begin;
                T inner_integral = T{0};
                auto qcoordNL = qcoordNL_begin;
                for(size_t qNL = 0, qshiftNL = qshiftNL_begin; qNL < elNL->qnodes_count(); ++qNL, ++JNL, ++qcoordNL, ++qshiftNL) {
                    const T weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL) * mesh::jacobian(*JNL);
                    inner_integral += weight * temperature_in_qnodes[qshiftNL];
                }
                inner_integral *= elL->weight(qL);
                for(const size_t comp : std::ranges::iota_view{0, 2})
                    integral[comp] += inner_integral * (*dNdL)[comp];
            }
            return integral;
        };

        f *= p1;
        const T p2 = T{1} - p1;
#pragma omp parallel for default(none) shared(f, p2, mesh_proxy) firstprivate(integrate_temperature_nonloc)
        for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node) {
            std::array<T, 2> integral = {};
            for(const I eL : mesh_proxy.nodes_elements_map(node)) {
                const size_t iL = mesh_proxy.global_to_local_numbering(eL, node);
                for(const I eNL : mesh_proxy.neighbors(eL))
                    integral += integrate_temperature_nonloc(eL, eNL, iL);
            }
            for(const size_t comp : std::ranges::iota_view{0, 2})
                f[2 * node + comp] += p2 * integral[comp];
        }
    }
}

}

#endif
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
    const T nu = parameters.poisson(),
            E  = parameters.young();
    const T factor = parameters.alpha * E / (T{2} * (T{1} - nu));
    const std::array<std::vector<T>, 2> gradient = mesh_proxy.template gradient(parameters.delta_temperature);
    using namespace metamath::function;
    const std::array<std::vector<T>, 2> temperature_eps = { factor * mesh_proxy.approx_in_quad(gradient[0]),
                                                            factor * mesh_proxy.approx_in_quad(gradient[1]) };

    const auto integrate_temperature_loc = [&mesh_proxy](const std::array<std::vector<T>, 2>& temperature_eps, const size_t e, const size_t i) {
        std::array<T, 2> integral = {};
        const auto& el  = mesh_proxy.mesh().element_2d(e);
              auto  dNd = mesh_proxy.dNdX(e, i);
        for(size_t q = 0, quad_shift = mesh_proxy.quad_shift(e); q < el->qnodes_count(); ++q, ++dNd, ++quad_shift)
            for(size_t comp = 0; comp < 2; ++comp)
                integral[comp] += el->weight(q) * temperature_eps[comp][quad_shift] * (*dNd)[comp];
        return integral;
    };

#pragma omp parallel for default(none) shared(f, temperature_eps, p1, mesh_proxy, integrate_temperature_loc)
    for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node) {
        std::array<T, 2> integral = {};
        for(const I e : mesh_proxy.nodes_elements_map(node)) {
            const size_t i = mesh_proxy.global_to_local_numbering(e).find(node)->second;
            integral += integrate_temperature_loc(temperature_eps, e, i);
        }
        f[2 * node + X] += p1 * integral[X];
        f[2 * node + Y] += p1 * integral[Y];
    }

    if (p1 < MAX_NONLOCAL_WEIGHT<T>) {
        const auto integrate_temperature_nonloc = [&mesh_proxy](const std::array<std::vector<T>, 2>& temperature_eps,
                                                                const size_t eL, const size_t eNL, const size_t iL,
                                                                const Influence_Function& influence_function) {
            std::array<T, 2> integral = {};
            const auto& elL            = mesh_proxy.mesh().element_2d(eL ),
                      & elNL           = mesh_proxy.mesh().element_2d(eNL);
                  auto  dNdL           = mesh_proxy.dNdX(eL, iL);
                  auto  qcoordL        = mesh_proxy.quad_coord(eL);
            const auto  qcoordNL_begin = mesh_proxy.quad_coord(eNL);
            const auto  JNL_begin      = mesh_proxy.jacobi_matrix(eNL);
            const size_t quad_shiftNL_begin = mesh_proxy.quad_shift(eNL);
            for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
                auto JNL = JNL_begin;
                auto qcoordNL = qcoordNL_begin;
                std::array<T, 2> inner_integral = {};
                for(size_t qNL = 0, quad_shiftNL = quad_shiftNL_begin; qNL < elNL->qnodes_count(); ++qNL, ++JNL, ++qcoordNL, ++quad_shiftNL) {
                    const T weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL) * mesh_proxy.jacobian(*JNL);
                    inner_integral[X] += weight * temperature_eps[X][quad_shiftNL];
                    inner_integral[Y] += weight * temperature_eps[Y][quad_shiftNL];
                }
                integral[X] += elL->weight(qL) * (*dNdL)[X] * inner_integral[X];
                integral[Y] += elL->weight(qL) * (*dNdL)[Y] * inner_integral[Y];
            }
            return integral;
        };

        const T p2 = T{1} - p1;
#pragma omp parallel for default(none) shared(f, temperature_eps, p2, mesh_proxy, integrate_temperature_nonloc) firstprivate(influence_function)
        for(size_t node = mesh_proxy.first_node(); node < mesh_proxy.last_node(); ++node) {
            std::array<T, 2> integral = {};
            for(const I eL : mesh_proxy.nodes_elements_map(node)) {
                const size_t iL = mesh_proxy.global_to_local_numbering(eL).find(node)->second;
                for(const I eNL : mesh_proxy.neighbors(eL))
                    integral += integrate_temperature_nonloc(temperature_eps, eL, eNL, iL, influence_function);
            }
            f[2 * node + X] += p2 * integral[X];
            f[2 * node + Y] += p2 * integral[Y];
        }
    }
}

}

#endif
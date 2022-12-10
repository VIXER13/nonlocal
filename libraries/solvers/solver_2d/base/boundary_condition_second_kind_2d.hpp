#ifndef NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP

#include "boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include "mesh_2d.hpp"

#include <eigen3/Eigen/Dense>

namespace nonlocal {

template<size_t DoF, class T, class I, class Conditions_Map>
void boundary_condition_second_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const mesh::mesh_2d<T, I>& mesh,
                                       const Conditions_Map& boundaries_conditions) {
    const auto integrate = [&mesh = mesh.container()](const auto& condition, const size_t be, const size_t i) {
        T integral = T{0};
        const auto el_data = mesh.element_1d_data(be);
        const auto& [_, __, el] = el_data;
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * condition(el_data.quad_coord(q)) * mesh::jacobian(el_data.jacobi_matrix(q));
        return integral;
    };

    utils::run_by_boundaries<DoF, second_kind_2d>(mesh.container(), boundaries_conditions,
        [&f, &integrate, &mesh, process_nodes = mesh.process_nodes()](const auto& condition, const size_t be, const size_t node, const size_t degree) {
            if (node >= process_nodes.front() && node <= process_nodes.back()) {
                const size_t index = DoF * node + degree;
                f[index] += integrate(condition, be, mesh.global_to_local(be, node));
            }
        });
}

}

#endif
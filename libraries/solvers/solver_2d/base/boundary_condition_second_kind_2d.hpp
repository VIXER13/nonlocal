#ifndef NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP

#include "boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include "mesh_2d.hpp"

#include <Eigen/Dense>

namespace nonlocal {

template<class T, class I, physics_t Physics, size_t DoF>
void boundary_condition_second_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const mesh::mesh_2d<T, I>& mesh,
                                       const boundaries_conditions_2d<T, Physics, DoF>& boundaries_conditions) {
    static constexpr auto integrate = [](const second_kind_2d<T, Physics>& condition, const auto& element, const size_t i) {
        T integral = T{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * condition(element.quad_coord(q)) * mesh::jacobian(element.jacobi_matrix(q));
        return integral;
    };

    utils::run_by_boundaries<second_kind_2d, Physics>(mesh.container(), boundaries_conditions,
        [&f, &mesh, process_nodes = mesh.process_nodes()]
        (const second_kind_2d<T, Physics>& condition, const size_t be, const size_t node, const size_t degree) {
            if (node >= process_nodes.front() && node <= process_nodes.back()) {
                const size_t index = DoF * node + degree;
                f[index] += integrate(condition, mesh.container().element_1d_data(be), mesh.global_to_local(be, node));
            }
        });
}

}

#endif
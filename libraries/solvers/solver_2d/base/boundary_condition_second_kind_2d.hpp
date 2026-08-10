#pragma once

#include "boundary_conditions_2d.hpp"
#include "solvers_utils.hpp"

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d {

template<class T, std::floating_point U, physics_t Physics, size_t DoF>
    void boundary_condition_second_kind_2d(std::vector<T>& f, const mesh::mesh_2d<U>& mesh,
                                           const boundaries_conditions_2d<U, Physics, DoF>& boundaries_conditions) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>>, "The floating point type of the mesh and the vector shall be the same");

    static constexpr auto integrate = [](const second_kind_2d<U, Physics>& condition, const auto& element, const size_t i) {
        U integral = U{0};
        const auto& [mesh, be] = element;
        const auto& el = mesh.element_1d(be);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()})
            integral += el.weight(q) * el.qN(i, q) * condition(element.quad_coord(q)) * mesh::jacobian(element.jacobi_matrix(q));
        return integral;
    };

    utils::run_by_boundaries<second_kind_2d, Physics>(mesh.container(), boundaries_conditions,
        [&f, &mesh, process_nodes = mesh.process_nodes()]
        (const second_kind_2d<U, Physics>& condition, const size_t be, const size_t node, const size_t degree) {
            if (node >= process_nodes.front() && node <= process_nodes.back()) {
                if constexpr (DoF == 1)
                    f[node] += integrate(condition, mesh.container().element_1d_data(be), mesh.global_to_local(be, node));
                else {
                    using namespace metamath::operators;
                    f[node][degree] += integrate(condition, mesh.container().element_1d_data(be), mesh.global_to_local(be, node));
                }
            }
        });
}

}
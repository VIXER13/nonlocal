#pragma once

#include <mesh/mesh_2d/mesh_2d.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace nonlocal {

template<size_t DoF, class T, class I, class Functor>
void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& right_part,
                          const mesh::mesh_2d<T, I>& mesh,
                          const Functor& functor) {
    const auto integrate = [&mesh, &functor](const size_t e, const size_t i) {
        std::conditional_t<DoF == 1, T, std::array<T, DoF>> integral = {};
        const auto& el = mesh.container().element_2d(e);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
            using namespace metamath::functions;
            integral += el.weight(q) * el.qN(i, q) * mesh::jacobian(mesh.jacobi_matrix(e, q)) * functor(mesh.quad_coord(e, q));
        }
        return integral;
    };

    const auto process_nodes = mesh.process_nodes();
#pragma parallel for default(none) shared(right_part, mesh, process_nodes, integrate)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        const size_t index = DoF * (node - process_nodes.front());
        for(const size_t e : mesh.elements(node)) {
            const size_t i = mesh.global_to_local(e, node);
            const std::array<T, DoF> integral = {integrate(e, i)};
            for(const size_t degree : std::ranges::iota_view{0u, DoF})
                right_part[index + degree] += integral[degree];
        }
    }
}

}
#pragma once

#include <mesh/mesh_2d/mesh_2d.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace nonlocal {

template<size_t DoF, class T, class I, class Functor>
void initialize_solution(Eigen::Matrix<T, Eigen::Dynamic, 1>& solution,
                         const mesh::mesh_2d<T, I>& mesh,
                         const Functor& functor) {
    const auto process_nodes = mesh.process_nodes();
#pragma parallel for default(none) shared(right_part, mesh, process_nodes, integrate)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node) {
        const size_t index = DoF * (node - process_nodes.front());
        const auto coord = mesh.container().node_coord(index);
        for(const size_t degree : std::ranges::iota_view{0u, DoF})
            solution[index + degree] = functor(coord);
    }
}

}
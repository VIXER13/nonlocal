#pragma once

#include <mesh/mesh_2d/mesh_2d.hpp>

namespace nonlocal::solver_2d {

template<class T, std::floating_point U, class Functor>
void integrate_right_part(std::vector<T>& right_part, const mesh::mesh_2d<U>& mesh, const Functor& functor) {
    static_assert(std::is_same_v<U, metamath::types::container_type_t<T>>, "The floating point type of the mesh and the vector shall be the same");

    const auto integrate = [&mesh, &functor](const size_t e, const size_t i) {
        T integral = {};
        const auto& el = mesh.container().element_2d(e);
        for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
            using namespace metamath::operators;
            integral += el.weight(q) * el.qN(i, q) * mesh.jacobian(e, q) * functor(mesh.quad_coord(e, q));
        }
        return integral;
    };

    using namespace metamath::operators;
    const auto process_nodes = mesh.process_nodes();
#pragma omp parallel for default(none) shared(right_part, mesh, process_nodes, integrate)
    for(size_t node = process_nodes.front(); node < *process_nodes.end(); ++node)
        for(const size_t e : mesh.elements(node))
            right_part[node - process_nodes.front()] += integrate(e, mesh.global_to_local(e, node));
}

}
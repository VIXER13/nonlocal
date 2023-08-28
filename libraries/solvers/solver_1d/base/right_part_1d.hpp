#ifndef NONLOCAL_RIGHT_PART_1D_HPP
#define NONLOCAL_RIGHT_PART_1D_HPP

#include "mesh_1d.hpp"

namespace nonlocal {

template<class T, class Vector, class Right_Part>
void integrate_right_part(Vector& f, const mesh::mesh_1d<T>& mesh, const Right_Part& right_part) {
    const auto integrate_function_on_element = [&mesh, &func = right_part](const size_t e, const size_t i) {
        T integral = T{0};
        const auto& el = mesh.element();
        for(const size_t q : std::ranges::iota_view{size_t{0}, el.qnodes_count()})
            integral += el.weight(q) * el.qN(0, i, q) * func(mesh.qnode_coord(e, q));
        return integral * mesh.jacobian(mesh.segment_number(e));
    };
#pragma omp parallel for default(none) shared(f, mesh, integrate_function_on_element)
    for(size_t node = 0; node < mesh.nodes_count(); ++node)
        for(const auto node_data : mesh.node_elements(node).to_array())
            if (node_data)
                f[node] -= integrate_function_on_element(node_data.element, node_data.node);
}

}

#endif
#ifndef RIGHT_PART_1D_HPP
#define RIGHT_PART_1D_HPP

#include "mesh_1d.hpp"
#include "boundary_condition_1d.hpp"
#include <array>
#include <vector>
#include <utility>
#include <unordered_map>

namespace nonlocal {

template<class T, class Vector>
void boundary_condition_first_kind_1d(Vector& f, const std::array<stationary_boundary_1d_t<T>, 2>& boundary_condition,
                                      const std::array<std::unordered_map<size_t, T>, 2>& matrix_bound) {
    const std::array<size_t, 2> ind = {0, size_t(f.size()-1)};
    for(size_t b = 0; b < boundary_condition.size(); ++b)
        if (boundary_type(boundary_condition[b]) == boundary_condition_t::FIRST_KIND) {
            for(const auto& [i, val] : matrix_bound[b])
                f[i] -= val * boundary_value(boundary_condition[b]);
            f[ind[b]] = boundary_value(boundary_condition[b]);
        }
}

template<class T, class Vector>
void boundary_condition_second_kind_1d(Vector& f, const std::array<stationary_boundary_1d_t<T>, 2>& boundary_condition,
                                       const std::array<size_t, 2>& ind) {
    for(size_t b = 0; b < boundary_condition.size(); ++b)
        if (boundary_type(boundary_condition[b]) == boundary_condition_t::SECOND_KIND ||
            boundary_type(boundary_condition[b]) == boundary_condition_t::THIRD_KIND)
            f[ind[b]] += boundary_value(boundary_condition[b]);
}

template<class T, class Function>
T integrate_function_on_element(const mesh::mesh_1d<T>& mesh, const size_t e, const size_t i, const Function& func) {
    T integral = T{0};
    const auto& el = mesh.element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * func(mesh.quad_coord(e, q));
    return integral * mesh.jacobian();
}

template<class T, class Vector, class Right_Part>
void integrate_right_part(Vector& f, const mesh::mesh_1d<T>& mesh, const Right_Part& right_part) {
#pragma omp parallel for default(none) shared(f, mesh, right_part)
    for(size_t node = 0; node < mesh.nodes_count(); ++node)
        for(const auto& [e, i] : mesh.node_elements(node).arr)
            if (e != std::numeric_limits<size_t>::max())
                f[node] += integrate_function_on_element(mesh, e, i, right_part);
}

}

#endif
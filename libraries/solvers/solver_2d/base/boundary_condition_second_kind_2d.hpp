#ifndef NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP
#define NONLOCAL_BOUNDARY_CONDITION_SECOND_KIND_2D_HPP

#include "boundary_conditions_2d.hpp"
#include "mesh_2d.hpp"

namespace nonlocal {

// template<class B, class T, class I, size_t DoF>
// void boundary_condition_second_kind_2d(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
//                                        const mesh::mesh_proxy<T, I>& mesh_proxy,
//                                        const std::unordered_map<std::string, stationary_boundary_2d_t<B, T, DoF>>& boundary_condition) {
//     const auto& mesh = mesh_proxy.mesh();
//     const std::array<size_t, 2>& first_last = {DoF * mesh_proxy.first_node(), DoF * mesh_proxy.last_node()};
//     for(const std::string& b : mesh.boundary_names())
//         for(const size_t comp : std::views::iota(size_t{0}, DoF))
//             if(boundary_condition.at(b)[comp].type == B(boundary_condition_t::SECOND_KIND))
//                 for(const size_t e : std::views::iota(size_t{0}, mesh.elements_count(b)))
//                     for(const size_t i : std::views::iota(size_t{0}, mesh.element_1d(b, e)->nodes_count()))
//                         if(const I node = DoF * mesh.node_number(b, e, i) + comp;
//                            node >= first_last.front() && node < first_last.back())
//                             f[node - first_last.front()] += integrate_boundary_gradient(mesh_proxy, b, e, i, boundary_condition.at(b)[comp].func);
// }

// template<class T, class I, class Boundary_Gradient>
// T integrate_boundary_gradient(const mesh::mesh_proxy<T, I>& mesh_proxy,
//                               const std::string& b, const size_t e, const size_t i,
//                               const Boundary_Gradient& boundary_gradient) {
//     T integral = 0;
//     const auto& be = mesh_proxy.mesh().element_1d(b, e);
//     auto qcoord = mesh_proxy.quad_coord(b, e);
//     auto J      = mesh_proxy.jacobi_matrix(b, e);
//     for(size_t q = 0; q < be->qnodes_count(); ++q, ++qcoord, ++J)
//         integral += be->weight(q) * be->qN(i, q) * boundary_gradient(*qcoord) * mesh::jacobian(*J);
//     return integral;
// }

// template<class T, class I, size_t DoF, class Vector>
// void boundary_condition_second_kind_2d(Vector& f,
//                                        const mesh::mesh_2d<T, I>& mesh,
//                                        const std::string& boundary_name,
//                                        const boundary_conditions_2d<T, DoF>& boundary_conditions) {
//     for(const size_t comp : std::ranges::iota_view{0u, DoF})
//         if(const auto* const condition = dynamic_cast<second_kind_2d<T>*>(boundary_conditions[comp].get()))
//             for(const size_t e : std::ranges::iota_view{0u, mesh.elements_count()})
//                 for(const size_t i : std::views::iota{0u, mesh.element_1d(b, e)->nodes_count()})
//                     if(const I node = DoF * mesh.node_number(b, e, i) + comp; node >= first_last.front() && node < first_last.back())
//                         f[node - first_last.front()] += 
// }

}

#endif
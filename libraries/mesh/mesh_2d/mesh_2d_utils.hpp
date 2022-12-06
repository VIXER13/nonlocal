#ifndef NONLOCAL_MESH_UTILS_HPP
#define NONLOCAL_MESH_UTILS_HPP

#include "mesh_2d.hpp"

namespace nonlocal::mesh::utils {

// template<class T, class I, std::ranges::random_access_range Vector>
// T integrate(const mesh_proxy<T, I>& mesh, const Vector& x) {
//     if (mesh.mesh().nodes_count() != size_t(x.size()))
//         throw std::logic_error{"The integral cannot be found because the vector size does not match the number of nodes."};
//     T integral = 0;
//     for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
//         auto J = mesh.jacobi_matrix(e);
//         const auto& el = mesh.mesh().element_2d(e);
//         for(size_t q = 0; q < el->qnodes_count(); ++q, ++J)
//             for(const size_t i : std::views::iota(size_t{0}, el->nodes_count()))
//                 integral += el->weight(q) * x[mesh.mesh().node_number(e, i)] * el->qN(i, q) * jacobian(*J);
//     }
//     return integral;
// }

template<class T, class I, class Vector>
std::array<std::vector<T>, 2> gradient_in_qnodes(const mesh_2d<T, I>& mesh, const Vector& x) {
    if (mesh.container().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    const size_t quadratures_count = mesh.quad_shift(mesh.container().elements_2d_count());
    std::array<std::vector<T>, 2> gradient{
        std::vector<T>(quadratures_count, T{0}),
        std::vector<T>(quadratures_count, T{0})
    };
#pragma parallel for default(none) shared(gradient, mesh, x)
    for(size_t e = 0; e < mesh.container().elements_2d_count(); ++e) {
        const auto& el = mesh.container().element_2d(e);
        for(size_t q = 0, qshift = mesh.quad_shift(e); q < el.qnodes_count(); ++q, ++qshift) {
            for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()}) {
                const std::array<T, 2>& derivatives = mesh.derivatives(e, i, q);
                const T& val = x[mesh.container().node_number(e, i)];
                gradient[X][qshift] += val * derivatives[X];
                gradient[Y][qshift] += val * derivatives[Y];
            }
            const T jac = jacobian(mesh.jacobi_matrix(qshift));
            gradient[X][qshift] /= jac;
            gradient[Y][qshift] /= jac;
        }
    }
    return gradient;
}

template<class T, class I, class Vector>
std::vector<T> qnodes_to_nodes(const mesh_2d<T, I>& mesh, const Vector& x) {
    if (mesh.quad_shift(mesh.container().elements_2d_count()) != size_t(x.size()))
        throw std::logic_error{"Cannot approximate node values because vector size does not match number of quadrature nodes"};
    std::vector<T> approximation(mesh.container().nodes_count(), T{0});
#pragma parallel for default(none) shared(approximation, mesh, x)
    for(const size_t node : mesh.container().nodes()) {
        T node_area = T{0};
        for(const I e : mesh.elements(node)) {
            const T area = mesh.element_area(e);
            const size_t i = mesh.global_to_local(e, node);
            const auto& el = mesh.container().element_2d(e);
            const size_t qshift = mesh.quad_shift(e) + el.nearest_qnode(i);
            approximation[node] += area * x[qshift];
            node_area += area;
        }
        approximation[node] /= node_area;
    }
    return approximation;
}

// template<class T, class I, std::ranges::random_access_range Vector>
// std::vector<T> approximate_in_qnodes(const mesh_proxy<T, I>& mesh, const Vector& x) {
//     if (mesh.mesh().nodes_count() != size_t(x.size()))
//         throw std::logic_error{"The approximation in qnodes cannot be found because the vector size does not match the number of nodes."};
//     std::vector<T> approximation(mesh.quad_shift(mesh.mesh().elements_count()), T{0});
//     for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
//         const auto& el = mesh.mesh().element_2d(e);
//         for(size_t q = 0, qshift = mesh.quad_shift(e); q < el->qnodes_count(); ++q, ++qshift)
//             for(size_t i = 0; i < el->nodes_count(); ++i)
//                 approximation[qshift] += x[mesh.mesh().node_number(e, i)] * el->qN(i, q);
//     }
//     return approximation;
// }

}

#endif
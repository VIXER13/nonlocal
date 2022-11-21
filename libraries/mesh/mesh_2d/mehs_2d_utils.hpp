#ifndef MESH_UTILS_HPP
#define MESH_UTILS_HPP

#include "mesh_proxy.hpp"

#include <array>

namespace nonlocal::mesh {

template<class T, class I, std::ranges::random_access_range Vector>
T integrate(const mesh_proxy<T, I>& mesh, const Vector& x) {
    if (mesh.mesh().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The integral cannot be found because the vector size does not match the number of nodes."};
    T integral = 0;
    for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
        auto J = mesh.jacobi_matrix(e);
        const auto& el = mesh.mesh().element_2d(e);
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J)
            for(const size_t i : std::views::iota(size_t{0}, el->nodes_count()))
                integral += el->weight(q) * x[mesh.mesh().node_number(e, i)] * el->qN(i, q) * jacobian(*J);
    }
    return integral;
}

template<class T, class I, std::ranges::random_access_range Vector>
std::array<std::vector<T>, 2> approximate_gradient_in_qnodes(const mesh_proxy<T, I>& mesh, const Vector& x) {
    if (mesh.mesh().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The gradient cannot be found because the vector size does not match the number of nodes."};
    const size_t quadratures_count = mesh.quad_shift(mesh.mesh().elements_count());
    std::array<std::vector<T>, 2> gradient{
        std::vector<T>(quadratures_count, T{0}),
        std::vector<T>(quadratures_count, T{0})
    };
    for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
        auto J = mesh.jacobi_matrix(e);
        const auto& el = mesh.mesh().element_2d(e);
        for(size_t q = 0, qshift = mesh.quad_shift(e); q < el->qnodes_count(); ++q, ++qshift, ++J) {
            for(const size_t i : std::views::iota(size_t{0}, el->nodes_count())) {
                const auto dNd = std::next(mesh.dNdX(e, i), q);
                gradient[0][qshift] += (*dNd)[0] * x[mesh.mesh().node_number(e, i)];
                gradient[1][qshift] += (*dNd)[1] * x[mesh.mesh().node_number(e, i)];
            }
            const T jac = jacobian(*J);
            gradient[0][qshift] /= jac;
            gradient[1][qshift] /= jac;
        }
    }
    return gradient;
}

template<class T, class I, std::ranges::random_access_range Vector>
std::vector<T> from_qnodes_to_nodes(const mesh_proxy<T, I>& mesh, const Vector& x) {
    if (mesh.quad_shift(mesh.mesh().elements_count()) != x.size())
        throw std::logic_error{"Cannot approximate node values because vector size does not match number of quadrature nodes"};
    std::vector<T> approximation(mesh.mesh().nodes_count(), T{0});
    for(size_t node = 0; node < mesh.mesh().nodes_count(); ++node) {
        T node_area = T{0};
        for(const I e : mesh.nodes_elements_map(node)) {
            const auto& el = mesh.mesh().element_2d(e);
            const size_t i = mesh.global_to_local_numbering(e, node);
            const size_t qshift = mesh.quad_shift(e) + el->nearest_qnode(i);
            const auto J = std::next(mesh.jacobi_matrix(e), el->nearest_qnode(i));
            approximation[node] += x[qshift] * mesh.element_area(e);
            node_area += mesh.element_area(e);
        }
        approximation[node] /= node_area;
    }
    return approximation;
}

template<class T, class I, std::ranges::random_access_range Vector>
std::vector<T> approximate_in_qnodes(const mesh_proxy<T, I>& mesh, const Vector& x) {
    if (mesh.mesh().nodes_count() != size_t(x.size()))
        throw std::logic_error{"The approximation in qnodes cannot be found because the vector size does not match the number of nodes."};
    std::vector<T> approximation(mesh.quad_shift(mesh.mesh().elements_count()), T{0});
    for(size_t e = 0; e < mesh.mesh().elements_count(); ++e) {
        const auto& el = mesh.mesh().element_2d(e);
        for(size_t q = 0, qshift = mesh.quad_shift(e); q < el->qnodes_count(); ++q, ++qshift)
            for(size_t i = 0; i < el->nodes_count(); ++i)
                approximation[qshift] += x[mesh.mesh().node_number(e, i)] * el->qN(i, q);
    }
    return approximation;
}

}

#endif
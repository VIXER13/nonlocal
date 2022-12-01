#ifndef NONLOCAL_MESH_CONTAINER_2D_UTILS_HPP
#define NONLOCAL_MESH_CONTAINER_2D_UTILS_HPP

#include "mesh_container_2d.hpp"

namespace nonlocal::mesh::utils {

template<class T, class I>
std::vector<std::vector<I>> node_elements_2d(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::vector<I>> node_elements(mesh.nodes_count());
    for(const size_t e : mesh.elements_2d()) {
        for(const size_t node : mesh.nodes(e))
            node_elements[node].push_back(e);
    }
    for(std::vector<I>& elements : node_elements)
        elements.shrink_to_fit();
    return node_elements;
}

template<class T, class I>
std::vector<std::unordered_map<I, uint8_t>> global_to_local(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh.elements_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
    for(size_t e = 0; e < mesh.elements_2d_count(); ++e) {
        const std::vector<I>& nodes = mesh.nodes(e);
        for(const size_t i : std::ranges::iota_view{0u, nodes.size()})
            global_to_local_numbering[e][nodes[i]] = i;
    }
    return global_to_local_numbering;
}

template<class T, class I, class Shift>
std::vector<I> quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh, const Shift& shift) {
    std::vector<I> quad_shifts(mesh.elements_2d_count() + 1);
    quad_shifts[0] = 0;
    for(const size_t e : mesh.elements_2d())
        quad_shifts[e + 1] = quad_shifts[e] + shift(e);
    return quad_shifts;
}

template<class T, class I>
std::vector<I> elements_quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh) {
    return quadrature_shifts_2d(mesh, [&mesh](const size_t e) { return mesh.element_2d(e).qnodes_count(); });
}

template<class T, class I>
std::vector<I> element_node_shits_quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh) {
    return quadrature_shifts_2d(mesh, [&mesh](const size_t e) { 
        const auto& el = mesh.element_2d(e);
        return el.nodes_count() * el.qnodes_count();
    });
}

template<template<class, size_t> class Output, class T, class I, class Functor>
std::vector<Output<T, 2>> approx_in_all_quad_nodes(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts, const Functor& functor) {
    if(mesh.elements_2d_count() + 1 != qshifts.size())
        throw std::logic_error{"The number of quadrature shifts and elements does not match."};
    std::vector<Output<T, 2>> data(qshifts.back());
    for(const size_t e : mesh.elements_2d()) {
        const auto element_data = mesh.element_2d_data(e);
        for(const size_t q : std::ranges::iota_view{0u, element_data.element.qnodes_count()})
            data[qshifts[e] + q] = functor(element_data, q);
    }
    return data;
}

template<class T, class I>
std::vector<std::array<T, 2>> approx_all_quad_nodes(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts) {
    return approx_in_all_quad_nodes<std::array>(mesh, qshifts, 
        [](const auto& element_data, const size_t q) { return element_data.quad_coord(q); });
}

template<class T, class I>
std::vector<metamath::types::square_matrix<T, 2>> approx_all_jacobi_matrices(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts) {
    return approx_in_all_quad_nodes<metamath::types::square_matrix>(mesh, qshifts, 
        [](const auto& element_data, const size_t q) { return element_data.jacobi_matrix(q); });
}

template<class T, class I>
std::vector<std::array<T, 2>> derivatives_in_quad(const mesh_container_2d<T, I>& mesh,
                                                  const std::vector<I>& quad_element_shifts,
                                                  const std::vector<I>& quad_nodes_shifts,
                                                  const std::vector<metamath::types::square_matrix<T, 2>>& jacobi_matrices) {
    if (mesh.elements_2d_count() + 1 != quad_element_shifts.size() || mesh.elements_2d_count() + 1 != quad_nodes_shifts.size())
        throw std::logic_error{"The number of quadrature shifts and elements does not match."};
    if (quad_element_shifts.back() != jacobi_matrices.size())
        throw std::logic_error{"The size of Jacobi matrices vector does not match with the quadratures nodes count."};
    std::vector<std::array<T, 2>> derivatives(quad_nodes_shifts.back());
#pragma omp parallel for default(none) shared(mesh, quad_element_shifts, quad_nodes_shifts, jacobi_matrices, derivatives)
    for(size_t e = 0; e < mesh.elements_2d_count(); ++e) {
        const auto& el = mesh.element_2d(e);
        for(const size_t i : std::ranges::iota_view{0u, el.nodes_count()})
            for(const size_t q : std::ranges::iota_view{0u, el.qnodes_count()}) {
                const metamath::types::square_matrix<T, 2>& J = jacobi_matrices[quad_element_shifts[e] + q];
                derivatives[quad_nodes_shifts[e] + i * el.qnodes_count() + q] = {
                     el.qNxi(i, q) * J[1][1] - el.qNeta(i, q) * J[1][0],
                    -el.qNxi(i, q) * J[0][1] + el.qNeta(i, q) * J[0][0]
                };
            }
    }
    return derivatives;
}

}

#endif
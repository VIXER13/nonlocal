#ifndef NONLOCAL_MESH_CONTAINER_2D_UTILS_HPP
#define NONLOCAL_MESH_CONTAINER_2D_UTILS_HPP

#include "mesh_container_2d.hpp"

namespace nonlocal::mesh::utils {

template<class T, class I>
std::vector<std::vector<I>> node_elements_2d(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::vector<I>> node_elements(mesh.nodes_count());
    for(const size_t element : mesh.elements_2d()) {
        for(const size_t node : mesh.nodes(element))
            node_elements[node].push_back(element);
    }
    for(std::vector<I>& elements : node_elements)
        elements.shrink_to_fit();
    return node_elements;
}

template<class T, class I>
std::vector<std::unordered_map<I, uint8_t>> global_to_local(const mesh_container_2d<T, I>& mesh) {
    std::vector<std::unordered_map<I, uint8_t>> global_to_local_numbering(mesh.elements_count());
#pragma omp parallel for default(none) shared(mesh, global_to_local_numbering)
    for(const size_t element : mesh.elements_2d()) {
        const std::vector<I>& nodes = mesh.nodes(element);
        for(const size_t local_index : std::ranges::iota_view{0u, nodes.size()})
            global_to_local_numbering[element][nodes[local_index]] = local_index;
    }
    return global_to_local_numbering;
}

template<class T, class I>
std::vector<I> quadrature_shifts_2d(const mesh_container_2d<T, I>& mesh) {
    std::vector<I> quad_shifts(mesh.elements_2d_count() + 1);
    quad_shifts[0] = 0;
    for(const size_t element : mesh.elements_2d())
        quad_shifts[element + 1] = quad_shifts[element] + mesh.element_2d(element).qnodes_count();
    return quad_shifts;
}

template<template<class, size_t> class Output, class T, class I, class Functor>
std::vector<Output<T, 2>> approx_in_all_quad_nodes(const mesh_container_2d<T, I>& mesh, const std::vector<I>& qshifts, const Functor& functor) {
    if(mesh.elements_2d_count() + 1 != qshifts.size())
        throw std::logic_error{"Quadrature shifts are incorrect size."};
    std::vector<Output<T, 2>> data(qshifts.back());
    for(const size_t element : mesh.elements_2d()) {
        const auto element_data = mesh.element_2d_data(element);
        for(const size_t q : std::ranges::iota_view{0u, element_data.element.qnodes_count()})
            data[qshifts[element] + q] = functor(element_data, q);
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

}

#endif
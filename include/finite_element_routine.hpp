#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

#include "mesh_2d.hpp"

template<class Type, class Index>
void mesh_run_loc(const mesh_2d<Type, Index> &mesh, const std::function<void(size_t, size_t, size_t)> &rule) {
    const finite_element::element_2d_integrate_base<Type> *e = nullptr;
#pragma omp parallel for default(none) shared(mesh) firstprivate(rule, e)
    for(size_t el = 0; el < mesh.elements_count(); ++el) {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t j = 0; j < e->nodes_count(); ++j)
                rule(i, j, el);
    }
}

template<class Type, class Index>
std::vector<Index> quadrature_shifts_init(const mesh_2d<Type, Index> &mesh) {
    std::vector<Index> shifts(mesh.elements_count()+1);
    shifts[0] = 0;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
        shifts[el+1] = shifts[el] + mesh.element_2d(mesh.element_type(el))->qnodes_count();
    return shifts;
}

template<class Type, class Index>
void approx_quad_nodes_coords(const mesh_2d<Type, Index> &mesh, const finite_element::element_2d_integrate_base<Type> *e,
                              const size_t el, matrix<Type> &coords) {
    coords.resize(e->qnodes_count(), 2);
    memset(coords.data(), 0, coords.size() * sizeof(Type));
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t component = 0; component < 2; ++component)
                coords(q, component) += mesh.coord(mesh.node(el, i), component) * e->qN(i, q);
}

template<class Type, class Index>
matrix<Type> approx_all_quad_nodes_coords(const mesh_2d<Type, Index> &mesh, const std::vector<Index> &shifts) {
    assert(mesh.elements_count()+1 == shifts.size());
    matrix<Type> coords(shifts.back(), 2, 0.0);
    const finite_element::element_2d_integrate_base<Type> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el) {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            for(size_t i = 0; i < e->nodes_count(); ++i)
                for(size_t component = 0; component < 2; ++component)
                    coords(shifts[el]+q, component) += mesh.coord(mesh.node(el, i), component) * e->qN(i, q);
    }
    return coords;
}

template<class Type, class Index>
void approx_quad_nodes_coord_bound(const mesh_2d<Type, Index> &mesh, const finite_element::element_1d_integrate_base<Type> *be,
                                   const size_t b, const size_t el, matrix<Type> &coords) {
    coords.resize(be->qnodes_count(), 2);
    memset(coords.data(), 0, coords.size() * sizeof(Type));
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        for(size_t i = 0; i < be->nodes_count(); ++i)
            for(size_t comp = 0; comp < 2; ++comp)
                coords(q, comp) += mesh.coord(mesh.boundary(b)(el, i), comp) * be->qN(i, q);
}

template<class Type, class Index>
void approx_jacobi_matrices(const mesh_2d<Type, Index> &mesh, const finite_element::element_2d_integrate_base<Type> *e,
                            const size_t el, matrix<Type> &jacobi_matrices) {
    jacobi_matrices.resize(e->qnodes_count(), 4);
    memset(jacobi_matrices.data(), 0, jacobi_matrices.size() * sizeof(Type));
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            jacobi_matrices(q, 0) += mesh.coord(mesh.node(el, i), 0) * e->qNxi (i, q);
            jacobi_matrices(q, 1) += mesh.coord(mesh.node(el, i), 0) * e->qNeta(i, q);
            jacobi_matrices(q, 2) += mesh.coord(mesh.node(el, i), 1) * e->qNxi (i, q);
            jacobi_matrices(q, 3) += mesh.coord(mesh.node(el, i), 1) * e->qNeta(i, q);
        }
}

template<class Type, class Index>
matrix<Type> approx_all_jacobi_matrices(const mesh_2d<Type, Index> &mesh, const std::vector<Index> &shifts) {
    assert(mesh.elements_count()+1 == shifts.size());
    matrix<Type> jacobi_matrices(shifts.back(), 4, 0.0);
    const finite_element::element_2d_integrate_base<Type> *e = nullptr;
    for(size_t el = 0; el < mesh.elements_count(); ++el) {
        e = mesh.element_2d(mesh.element_type(el));
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            for(size_t i = 0; i < e->nodes_count(); ++i) {
                jacobi_matrices(shifts[el]+q, 0) += mesh.coord(mesh.node(el, i), 0) * e->qNxi (i, q);
                jacobi_matrices(shifts[el]+q, 1) += mesh.coord(mesh.node(el, i), 0) * e->qNeta(i, q);
                jacobi_matrices(shifts[el]+q, 2) += mesh.coord(mesh.node(el, i), 1) * e->qNxi (i, q);
                jacobi_matrices(shifts[el]+q, 3) += mesh.coord(mesh.node(el, i), 1) * e->qNeta(i, q);
            }
    }
    return jacobi_matrices;
}

template<class Type, class Index>
void approx_jacobi_matrices_bound(const mesh_2d<Type, Index> &mesh, const finite_element::element_1d_integrate_base<Type> *be,
                                  const size_t b, const size_t el, matrix<Type> &jacobi_matrices) {
    jacobi_matrices.resize(be->qnodes_count(), 2);
    memset(jacobi_matrices.data(), 0, jacobi_matrices.size() * sizeof(Type));
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        for(size_t i = 0; i < be->nodes_count(); ++i)
            for(size_t comp = 0; comp < 2; ++comp)
                jacobi_matrices(q, comp) += mesh.coord(mesh.boundary(b)(el, i), comp) * be->qNxi(i, q);
}

#endif
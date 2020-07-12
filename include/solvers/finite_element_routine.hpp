#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

#include <string.h>
#include <exception>
#include "mesh_2d.hpp"

// Rule - функтор с сигнатурой void(size_t, size_t, size_t)
template<class Type, class Index, class Rule>
static void mesh_run_loc(const mesh_2d<Type, Index>& mesh, const Rule& rule)
{
#pragma omp parallel for default(none) shared(mesh) firstprivate(rule)
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t j = 0; j < e->nodes_count(); ++j)
                rule(i, j, el);
    }
}

// Rule - функтор с сигнатурой void(size_t, size_t, size_t, size_t)
template<class Type, class Index, class Rule>
static void mesh_run_nonloc(const mesh_2d<Type, Index>& mesh, const Rule& rule)
{
#pragma omp parallel for default(none) shared(mesh) firstprivate(rule)
    for(size_t elL = 0; elL < mesh.elements_count(); ++elL)
    {
        const auto& eL = mesh.element_2d(mesh.element_type(elL));
        for(const auto elNL : mesh.neighbor(elL))
        {
            const auto& eNL = mesh.element_2d(mesh.element_type(elNL));
            for(size_t iL = 0; iL < eL->nodes_count(); ++iL)
                for(size_t jNL = 0; jNL < eNL->nodes_count(); ++jNL)
                    rule(iL, jNL, elL, elNL);
        }
    }
}

template<class Type, class Index>
static std::vector<Index> quadrature_shifts_init(const mesh_2d<Type, Index>& mesh)
{
    std::vector<Index> shifts(mesh.elements_count()+1);
    shifts[0] = 0;
    for(size_t el = 0; el < mesh.elements_count(); ++el)
        shifts[el+1] = shifts[el] + mesh.element_2d(mesh.element_type(el))->qnodes_count();
    return shifts;
}

template<class Type, class Index, class Finite_Element_2D_Pointer>
static void approx_quad_nodes_coords(const mesh_2d<Type, Index>& mesh, const Finite_Element_2D_Pointer& e,
                                     const size_t el, matrix<Type>& coords)
{
    coords.resize(e->qnodes_count(), 2);
    memset(coords.data(), 0, coords.size() * sizeof(Type));
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
            for(size_t component = 0; component < 2; ++component)
                coords(q, component) += mesh.coord(mesh.node_number(el, i), component) * e->qN(i, q);
}

template<class Type, class Index>
static matrix<Type> approx_all_quad_nodes_coords(const mesh_2d<Type, Index>& mesh, const std::vector<Index>& shifts)
{
    if(mesh.elements_count()+1 != shifts.size())
        throw std::logic_error("mesh.elements_count()+1 != shifts.size()");
    matrix<Type> coords(shifts.back(), 2, 0.);
#pragma omp parallel for default(none) shared(mesh, shifts, coords)
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            for(size_t i = 0; i < e->nodes_count(); ++i)
                for(size_t component = 0; component < 2; ++component)
                    coords(shifts[el]+q, component) += mesh.coord(mesh.node_number(el, i), component) * e->qN(i, q);
    }
    return coords;
}

template<class Type, class Index, class Finite_Element_1D_Pointer>
static void approx_quad_nodes_coord_bound(const mesh_2d<Type, Index>& mesh, const Finite_Element_1D_Pointer& be,
                                          const size_t b, const size_t el, matrix<Type>& coords)
{
    coords.resize(be->qnodes_count(), 2);
    memset(coords.data(), 0, coords.size() * sizeof(Type));
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        for(size_t i = 0; i < be->nodes_count(); ++i)
            for(size_t comp = 0; comp < 2; ++comp)
                coords(q, comp) += mesh.coord(mesh.boundary(b)(el, i), comp) * be->qN(i, q);
}

template<class Type, class Index, class Finite_Element_2D_Pointer>
static void approx_jacobi_matrices(const mesh_2d<Type, Index>& mesh, const Finite_Element_2D_Pointer& e,
                                   const size_t el, matrix<Type>& jacobi_matrices) 
{
    jacobi_matrices.resize(e->qnodes_count(), 4);
    memset(jacobi_matrices.data(), 0, jacobi_matrices.size() * sizeof(Type));
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        for(size_t i = 0; i < e->nodes_count(); ++i)
        {
            jacobi_matrices(q, 0) += mesh.coord(mesh.node_number(el, i), 0) * e->qNxi (i, q);
            jacobi_matrices(q, 1) += mesh.coord(mesh.node_number(el, i), 0) * e->qNeta(i, q);
            jacobi_matrices(q, 2) += mesh.coord(mesh.node_number(el, i), 1) * e->qNxi (i, q);
            jacobi_matrices(q, 3) += mesh.coord(mesh.node_number(el, i), 1) * e->qNeta(i, q);
        }
}

template<class Type, class Index>
static matrix<Type> approx_all_jacobi_matrices(const mesh_2d<Type, Index>& mesh, const std::vector<Index>& shifts)
{
    if(mesh.elements_count()+1 != shifts.size())
        throw std::logic_error("mesh.elements_count()+1 != shifts.size()");
    matrix<Type> jacobi_matrices(shifts.back(), 4, 0);
#pragma omp parallel for default(none) shared(mesh, shifts, jacobi_matrices)
    for(size_t el = 0; el < mesh.elements_count(); ++el)
    {
        const auto& e = mesh.element_2d(mesh.element_type(el));
        for(size_t q = 0; q < e->qnodes_count(); ++q)
            for(size_t i = 0; i < e->nodes_count(); ++i)
            {
                jacobi_matrices(shifts[el]+q, 0) += mesh.coord(mesh.node_number(el, i), 0) * e->qNxi (i, q);
                jacobi_matrices(shifts[el]+q, 1) += mesh.coord(mesh.node_number(el, i), 0) * e->qNeta(i, q);
                jacobi_matrices(shifts[el]+q, 2) += mesh.coord(mesh.node_number(el, i), 1) * e->qNxi (i, q);
                jacobi_matrices(shifts[el]+q, 3) += mesh.coord(mesh.node_number(el, i), 1) * e->qNeta(i, q);
            }
    }
    return jacobi_matrices;
}

template<class Type, class Index, class Finite_Element_1D_Pointer>
static void approx_jacobi_matrices_bound(const mesh_2d<Type, Index>& mesh, const Finite_Element_1D_Pointer& be,
                                         const size_t b, const size_t el, matrix<Type>& jacobi_matrices)
{
    jacobi_matrices.resize(be->qnodes_count(), 2);
    memset(jacobi_matrices.data(), 0, jacobi_matrices.size() * sizeof(Type));
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        for(size_t i = 0; i < be->nodes_count(); ++i)
            for(size_t comp = 0; comp < 2; ++comp)
                jacobi_matrices(q, comp) += mesh.coord(mesh.boundary(b)(el, i), comp) * be->qNxi(i, q);
}

// Right_Part - функтор с сигнатурой Type(Type, Type)
template<class Type, class Finite_Element_2D_Pointer, class Right_Part>
static Type integrate_right_part_function(const Finite_Element_2D_Pointer& e, const size_t i,
                                          const matrix<Type>& coords, const matrix<Type>& jacobi_matrices,
                                          const Right_Part& fun)
{
    Type integral = 0;
    for(size_t q = 0; q < e->qnodes_count(); ++q)
        integral += e->weight(q) * e->qN(i, q) * fun(coords(q, 0), coords(q, 1)) *
                    (jacobi_matrices(q, 0)*jacobi_matrices(q, 3) - jacobi_matrices(q, 1)*jacobi_matrices(q, 2));
    return integral;
}

// Boundary_Gradient - функтор с сигнатурой Type(Type, Type)
template<class Type, class Finite_Element_1D_Pointer, class Boundary_Gradient>
static Type integrate_boundary_gradient(const Finite_Element_1D_Pointer& be, const size_t i,
                                        const matrix<Type>& coords, const matrix<Type>& jacobi_matrices, 
                                        const Boundary_Gradient& boundary_gradient)
{
    Type integral = 0;
    for(size_t q = 0; q < be->qnodes_count(); ++q)
        integral += be->weight(q) * be->qN(i, q) * boundary_gradient(coords(q, 0), coords(q, 1)) *
                    sqrt(jacobi_matrices(q, 0)*jacobi_matrices(q, 0) + jacobi_matrices(q, 1)*jacobi_matrices(q, 1));
    return integral;
}

#endif
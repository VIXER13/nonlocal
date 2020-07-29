#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

// Базовые операции, которые требуются во всех конечно-элементных решателях.

#include "mesh/mesh.hpp"

namespace nonlocal {

enum class boundary_type : uint8_t {
    FIRST_KIND,
    SECOND_KIND
};

class _finite_element_routine {
protected:
    enum component : bool { X, Y };
    static constexpr double MAX_LOCAL_WEIGHT = 0.999;

    explicit _finite_element_routine() noexcept = default;

    // Функция обхода сетки в локальных постановках.
    // Нужна для предварительного подсчёта количества элементов и интегрирования системы.
    // Rule - функтор с сигнатурой void(size_t, size_t, size_t)
    template<class Type, class Index, class Rule>
    static void mesh_run_loc(const mesh::mesh_2d<Type, Index>& mesh, const Rule& rule) {
#pragma omp parallel for default(none) shared(mesh) firstprivate(rule)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t i = 0; i < e->nodes_count(); ++i)     // Проекционные функции
                for(size_t j = 0; j < e->nodes_count(); ++j) // Функции формы
                    rule(i, j, el);
        }
    }

    // Функция обхода сетки в нелокальных постановках.
    // Нужна для предварительного подсчёта количества элементов и интегрирования системы.
    // Rule - функтор с сигнатурой void(size_t, size_t, size_t, size_t)
    template<class Type, class Index, class Rule>
    static void mesh_run_nonloc(const mesh::mesh_2d<Type, Index>& mesh, const Rule& rule) {
#pragma omp parallel for default(none) shared(mesh) firstprivate(rule)
        for(size_t elL = 0; elL < mesh.elements_count(); ++elL) {
            const auto& eL = mesh.element_2d(mesh.element_2d_type(elL));
            for(const auto elNL : mesh.element_neighbors(elL)) {
                const auto& eNL = mesh.element_2d(mesh.element_2d_type(elNL));
                for(size_t iL = 0; iL < eL->nodes_count(); ++iL)         // Проекционные функции
                    for(size_t jNL = 0; jNL < eNL->nodes_count(); ++jNL) // Функции формы
                        rule(iL, jNL, elL, elNL);
            }
        }
    }

    template<class Type, class Index, class Rule>
    static void boundary_nodes_run(const mesh::mesh_2d<Type, Index>& mesh, 
                                   const Rule& rule) {
        for(size_t b = 0; b < mesh.boundary_groups_count(); ++b)
            for(size_t el = 0; el < mesh.elements_count(b); ++el) {
                const auto& be = mesh.element_1d(mesh.element_1d_type(b, el));
                for(size_t i = 0; i < be->nodes_count(); ++i)
                    rule(b, el, i);
            }
    }

    // Квадратурные сдвиги по элементам.
    template<class Type, class Index>
    static std::vector<Index> quadrature_shifts_init(const mesh::mesh_2d<Type, Index>& mesh) {
        std::vector<Index> shifts(mesh.elements_count()+1);
        shifts[0] = 0;
        for(size_t el = 0; el < mesh.elements_count(); ++el)
            shifts[el+1] = shifts[el] + mesh.element_2d(mesh.element_2d_type(el))->qnodes_count();
        return std::move(shifts);
    }

    // Аппроксимация глобальных координат всех квадратурных узлов сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    template<class Type, class Index>
    static std::vector<std::array<Type, 2>> approx_all_quad_nodes(const mesh::mesh_2d<Type, Index>& mesh,
                                                                  const std::vector<Index>& shifts) {
        if(mesh.elements_count()+1 != shifts.size())
            throw std::logic_error{"mesh.elements_count()+1 != shifts.size()"};
        std::vector<std::array<Type, 2>> coords(shifts.back(), std::array<Type, 2>{});
#pragma omp parallel for default(none) shared(mesh, shifts, coords)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                for(size_t i = 0; i < e->nodes_count(); ++i)
                    for(size_t comp = 0; comp < 2; ++comp)
                        coords[shifts[el]+q][comp] += mesh.node(mesh.node_number(el, i))[comp] * e->qN(i, q);
        }
        return std::move(coords);
    }

    template<class Type, class Index>
    static void approx_quad_nodes_on_bound(std::vector<std::array<Type, 2>>& quad_nodes, 
                                           const mesh::mesh_2d<Type, Index>& mesh,
                                           const size_t b, const size_t el) {
        const auto& be = mesh.element_1d(mesh.element_1d_type(b, el));
        quad_nodes.clear();
        quad_nodes.resize(be->qnodes_count(), {});
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            for(size_t i = 0; i < be->nodes_count(); ++i)
                for(size_t comp = 0; comp < 2; ++comp)
                    quad_nodes[q][comp] += mesh.node(mesh.node_number(b, el, i))[comp] * be->qN(i, q);
    }

    // Аппроксимация матриц Якоби во всех квадратурных узлах сетки.
    // Перед вызовом обязательно должны быть проинициализированы квадратурные сдвиги
    template<class Type, class Index>
    static std::vector<std::array<Type, 4>> approx_all_jacobi_matrices(const mesh::mesh_2d<Type, Index>& mesh, 
                                                                       const std::vector<Index>& shifts) {
        if(mesh.elements_count()+1 != shifts.size())
            throw std::logic_error{"mesh.elements_count()+1 != shifts.size()"};
        std::vector<std::array<Type, 4>> jacobi_matrices(shifts.back(), std::array<Type, 4>{});
#pragma omp parallel for default(none) shared(mesh, shifts, jacobi_matrices)
        for(size_t el = 0; el < mesh.elements_count(); ++el) {
            const auto& e = mesh.element_2d(mesh.element_2d_type(el));
            for(size_t q = 0; q < e->qnodes_count(); ++q)
                for(size_t i = 0; i < e->nodes_count(); ++i) {
                    jacobi_matrices[shifts[el]+q][0] += mesh.node(mesh.node_number(el, i))[0] * e->qNxi (i, q);
                    jacobi_matrices[shifts[el]+q][1] += mesh.node(mesh.node_number(el, i))[0] * e->qNeta(i, q);
                    jacobi_matrices[shifts[el]+q][2] += mesh.node(mesh.node_number(el, i))[1] * e->qNxi (i, q);
                    jacobi_matrices[shifts[el]+q][3] += mesh.node(mesh.node_number(el, i))[1] * e->qNeta(i, q);
                }
        }
        return std::move(jacobi_matrices);
    }

    template<class Type, class Index>
    static void approx_jacobi_matrices_on_bound(std::vector<std::array<Type, 2>>& jacobi_matrices, 
                                                const mesh::mesh_2d<Type, Index>& mesh,
                                                const size_t b, const size_t el) {
        const auto& be = mesh.element_1d(mesh.element_1d_type(b, el));
        jacobi_matrices.clear();
        jacobi_matrices.resize(be->qnodes_count(), {});
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            for(size_t i = 0; i < be->nodes_count(); ++i)
                for(size_t comp = 0; comp < 2; ++comp)
                    jacobi_matrices[q][comp] += mesh.node(mesh.node_number(b, el, i))[comp] * be->qNxi(i, q);
    }

    template<class T>
    static T jacobian(const std::array<T, 4>& jacobi_matrix) noexcept {
        return jacobi_matrix[0] * jacobi_matrix[3] - jacobi_matrix[1] * jacobi_matrix[2];
    }

    template<class T>
    static T jacobian(const std::array<T, 2>& jacobi_matrix) noexcept {
        return sqrt(jacobi_matrix[0] * jacobi_matrix[0] + jacobi_matrix[1] * jacobi_matrix[1]);
    }

    template<bool Component, class T, class Finite_Element_2D_Ptr>
    static T dNd(const Finite_Element_2D_Ptr& e, const size_t i, const size_t q,
                 const std::array<T, 4>& jacobi_matrix) {
        if constexpr (Component == component::X)
            return e->qNxi(i, q) * jacobi_matrix[3] - e->qNeta(i, q) * jacobi_matrix[2];
        if constexpr (Component == component::Y)
            return -e->qNxi(i, q) * jacobi_matrix[1] + e->qNeta(i, q) * jacobi_matrix[0];
    }

    // Function - функтор с сигнатурой T(std::array<T, 2>&)
    template<class T, class Finite_Element_2D_Ptr, class Function>
    static T integrate_function(const Finite_Element_2D_Ptr& e, const size_t i,
                                const std::vector<std::array<T, 2>>& quad_nodes,
                                const std::vector<std::array<T, 4>>& jacobi_matrices,
                                size_t quad_shift, const Function& func) {
        T integral = 0;
        for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
            integral += e->weight(q) * e->qN(i, q) * func(quad_nodes[quad_shift]) * jacobian(jacobi_matrices[quad_shift]);
        return integral;
    }

    // Boundary_Gradient - функтор с сигнатурой T(std::array<T, 2>&)
    template<class T, class Finite_Element_1D_Pointer, class Boundary_Gradient>
    static T integrate_boundary_gradient(const Finite_Element_1D_Pointer& be, const size_t i,
                                         const std::vector<std::array<T, 2>>& quad_nodes, 
                                         const std::vector<std::array<T, 2>>& jacobi_matrices, 
                                         const Boundary_Gradient& boundary_gradient) {
        T integral = 0;
        for(size_t q = 0; q < be->qnodes_count(); ++q)
            integral += be->weight(q) * be->qN(i, q) * boundary_gradient(quad_nodes[q]) * jacobian(jacobi_matrices[q]);
        return integral;
    }
};

}

#endif
#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

#include <numeric>
#include <petsc.h>
#include <petscsystypes.h>
#include "../../Eigen/Eigen/Sparse"
#undef I // for new version GCC, when use I macros
#include "mesh.hpp"
#include "right_partition.hpp"
#include "boundary_condition.hpp"

namespace nonlocal {

template<class T, class I>
class finite_element_solver_base {
    std::shared_ptr<mesh::mesh_info<T, I>> _mesh_info;

protected:
    using Finite_Element_1D_Ptr = typename mesh::mesh_2d<T, I>::Finite_Element_1D_Ptr;
    using Finite_Element_2D_Ptr = typename mesh::mesh_2d<T, I>::Finite_Element_2D_Ptr;

    enum component : bool {X, Y};
    static constexpr T MAX_LOCAL_WEIGHT = 0.999;

    explicit finite_element_solver_base(const std::shared_ptr<mesh::mesh_info<T, I>>& mesh) { set_mesh(mesh); }
    virtual ~finite_element_solver_base() noexcept = default;

    const mesh::mesh_2d<T, I>&            mesh                     ()                     const { return _mesh_info->mesh(); }
    int                                   rank                     ()                     const { return _mesh_info->rank(); }
    int                                   size                     ()                     const { return _mesh_info->size(); }
    size_t                                first_node               ()                     const { return _mesh_info->first_node(); }
    size_t                                last_node                ()                     const { return _mesh_info->last_node(); }
    I                                     quad_shift               (const size_t element) const { return _mesh_info->quad_shift(element); }
    const std::array<T, 2>&               quad_coord               (const size_t quad)    const { return _mesh_info->quad_coord(quad); }
    const std::array<T, 4>&               jacobi_matrix            (const size_t quad)    const { return _mesh_info->jacobi_matrix(quad); }
    const std::array<T, 4>&               jacobi_matrix_node       (const size_t quad)    const { return _mesh_info->jacobi_matrix_node(quad); }
    const std::vector<I>&                 nodes_elements_map       (const size_t node)    const { return _mesh_info->nodes_elements_map(node); }
    const std::unordered_map<I, uint8_t>& global_to_local_numbering(const size_t element) const { return _mesh_info->global_to_local_numbering(element); }
    const std::vector<I>&                 neighbors                (const size_t element) const { return _mesh_info->neighbors(element); }
    I                                     quad_shift               (const size_t bound, const size_t element) const { return _mesh_info->quad_shift(bound, element); }
    const std::array<T, 2>&               quad_coord               (const size_t bound, const size_t quad)    const { return _mesh_info->quad_coord(bound, quad); }
    const std::array<T, 2>&               jacobi_matrix            (const size_t bound, const size_t quad)    const { return _mesh_info->jacobi_matrix(bound, quad); }

    static T jacobian(const std::array<T, 4>& J) noexcept { return mesh::mesh_info<T, I>::jacobian(J);   }
           T jacobian(const size_t quad_shift)   const    { return _mesh_info->jacobian(quad_shift);     }
    static T jacobian(const std::array<T, 2>& J) noexcept { return std::sqrt(J[0] * J[0] + J[1] * J[1]); }

    template<bool Component>
    T dNd(const Finite_Element_2D_Ptr& e, const size_t i, const size_t q, const size_t quad_shift) const;

    // Функция обхода сетки в локальных постановках.
    // Нужна для предварительного подсчёта количества элементов и интегрирования системы.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t)
    template<class Callback>
    void mesh_run_loc(const Callback& callback) const;

    // Функция обхода сетки в нелокальных постановках.
    // Нужна для предварительного подсчёта количества элементов и интегрирования системы.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t, size_t)
    template<class Callback>
    void mesh_run_nonloc(const Callback& callback) const;

    // Функция обхода групп граничных элементов.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t)
    template<class Callback>
    void boundary_nodes_run(const Callback& callback) const;

    void convert_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                          const std::vector<std::unordered_set<I>>& portrait) const;

    // Function - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Function>
    T integrate_function(const Finite_Element_2D_Ptr& e,
                         const size_t i, size_t quad_shift,
                         const Function& func) const;

    template<size_t DoF>
    void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                              const right_partition<T, DoF>& right_part) const;

    // Boundary_Gradient - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Boundary_Gradient>
    T integrate_boundary_gradient(const Finite_Element_1D_Ptr& be,
                                  const size_t b, const size_t i, size_t quad_shift,
                                  const Boundary_Gradient& boundary_gradient) const;

    template<class B, size_t DoF>
    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const std::vector<boundary_condition<T, B, DoF>>& bounds_cond,
                                       const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound) const;

    template<class B, size_t DoF>
    void integrate_boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                  const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const;

    void PETSc_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                      const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K);

public:
    void set_mesh(const std::shared_ptr<mesh::mesh_info<T, I>>& mesh_info) { _mesh_info = mesh_info; }
    const std::shared_ptr<mesh::mesh_info<T, I>>& get_mesh_info() const { return _mesh_info; }
};

template<class T, class I>
template<bool Component>
T finite_element_solver_base<T, I>::dNd(const Finite_Element_2D_Ptr& e, const size_t i, const size_t q, const size_t quad_shift) const {
    const std::array<T, 4>& J = jacobi_matrix(quad_shift);
    if constexpr (Component == component::X)
        return  e->qNxi(i, q) * J[3] - e->qNeta(i, q) * J[2];
    if constexpr (Component == component::Y)
        return -e->qNxi(i, q) * J[1] + e->qNeta(i, q) * J[0];
}

template<class T, class I>
template<class Callback>
void finite_element_solver_base<T, I>::mesh_run_loc(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback)
    for(size_t node = first_node(); node < last_node(); ++node) {
        for(const I el : nodes_elements_map(node)) {
            const size_t i = global_to_local_numbering(el).find(node)->second; // Проекционные функции
            for(size_t j = 0; j < mesh().nodes_count(el); ++j) // Аппроксимационные функции
                callback(el, i, j);
        }
    }
}

template<class T, class I>
template<class Callback>
void finite_element_solver_base<T, I>::mesh_run_nonloc(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback) schedule(dynamic)
    for(size_t node = first_node(); node < last_node(); ++node) {
        for(const I elL : nodes_elements_map(node)) {
            const size_t iL = global_to_local_numbering(elL).find(node)->second; // Проекционные функции
            for(const I elNL : neighbors(elL)) {
                for(size_t jNL = 0; jNL < mesh().nodes_count(elNL); ++jNL) // Аппроксимационные функции
                    callback(elL, iL, elNL, jNL);
            }
        }
    }
}

template<class T, class I>
template<class Callback>
void finite_element_solver_base<T, I>::boundary_nodes_run(const Callback& callback) const {
    for(size_t b = 0; b < mesh().boundary_groups_count(); ++b)
        for(size_t el = 0; el < mesh().elements_count(b); ++el) {
            const auto& be = mesh().element_1d(b, el);
            for(size_t i = 0; i < be->nodes_count(); ++i)
                callback(b, el, i);
        }
}

template<class T, class I>
void finite_element_solver_base<T, I>::convert_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                                        const std::vector<std::unordered_set<I>>& portrait) const {
    static constexpr auto accumulator = [](const size_t sum, const std::unordered_set<I>& row) { return sum + row.size(); };
    K.data().resize(std::accumulate(portrait.cbegin(), portrait.cend(), size_t{0}, accumulator));
    K.outerIndexPtr()[0] = 0;
    for(size_t row = 0; row < portrait.size(); ++row)
        K.outerIndexPtr()[row+1] = K.outerIndexPtr()[row] + portrait[row].size();
#pragma omp parallel for default(none) shared(K, portrait)
    for(size_t row = 0; row < portrait.size(); ++row) {
        I inner_index = K.outerIndexPtr()[row];
        for(const I col : portrait[row]) {
            K.valuePtr()[inner_index] = 0;
            K.innerIndexPtr()[inner_index++] = col;
        }
        std::sort(&K.innerIndexPtr()[K.outerIndexPtr()[row]], &K.innerIndexPtr()[K.outerIndexPtr()[row+1]]);
    }
}

template<class T, class I>
template<class Function>
T finite_element_solver_base<T, I>::integrate_function(const Finite_Element_2D_Ptr& e,
                                                       const size_t i, size_t quad_shift,
                                                       const Function& func) const {
    T integral = 0;
    for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
        integral += e->weight(q) * e->qN(i, q) * func(quad_coord(quad_shift)) * jacobian(quad_shift);
    return integral;
}

template<class T, class I>
template<size_t DoF>
void finite_element_solver_base<T, I>::integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                            const right_partition<T, DoF>& right_part) const {
#pragma omp parallel for default(none) shared(f, right_part)
    for(size_t node = first_node(); node < last_node(); ++node)
        for(const I el : nodes_elements_map(node)) {
            const size_t i = global_to_local_numbering(el).find(node)->second;
            for(size_t comp = 0; comp < DoF; ++comp)
                f[DoF * (node - first_node()) + comp] += integrate_function(mesh().element_2d(el), i, quad_shift(el), right_part[comp]);
        }
}

template<class T, class I>
template<class Boundary_Gradient>
T finite_element_solver_base<T, I>::integrate_boundary_gradient(const Finite_Element_1D_Ptr& be,
                                                                const size_t b, const size_t i, size_t quad_shift,
                                                                const Boundary_Gradient& boundary_gradient) const {
    T integral = 0;
    for(size_t q = 0; q < be->qnodes_count(); ++q, ++quad_shift)
        integral += be->weight(q) * be->qN(i, q) * boundary_gradient(quad_coord(b, quad_shift)) * jacobian(jacobi_matrix(b, quad_shift));
    return integral;
}

template<class T, class I>
template<class B, size_t DoF>
void finite_element_solver_base<T, I>::boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                     const std::vector<boundary_condition<T, B, DoF>>& bounds_cond,
                                                                     const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(K_bound.cols());
    boundary_nodes_run(
        [this, &bounds_cond, &x](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < DoF; ++comp)
                if (bounds_cond[b].type(comp) == B(boundary_type::FIRST_KIND)) {
                    const I node = DoF * mesh().node_number(b, el, i) + comp;
                    if (x[node] == 0)
                        x[node] = bounds_cond[b].func(comp)(mesh().node(node));
                }
        });

    f -= K_bound * x;

    boundary_nodes_run(
        [this, &bounds_cond, &x, &f](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < DoF; ++comp)
                if (bounds_cond[b].type(comp) == B(boundary_type::FIRST_KIND)) {
                    const I node = DoF * mesh().node_number(b, el, i) + comp;
                    if (node >= DoF * first_node() && node < DoF * last_node())
                        f[node - DoF * first_node()] = x[node];
                }
        });
}

template<class T, class I>
template<class B, size_t DoF>
void finite_element_solver_base<T, I>::integrate_boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                                const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const {
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        for(size_t comp = 0; comp < DoF; ++comp)
            if(bounds_cond[b].type(comp) == B(boundary_type::SECOND_KIND)) {
                for(size_t el = 0; el < mesh().elements_count(b); ++el) {
                    const auto& be = mesh().element_1d(b, el);
                    for(size_t i = 0; i < be->nodes_count(); ++i) {
                        const I node = DoF * mesh().node_number(b, el, i) + comp;
                        if(node >= DoF * first_node() && node < DoF * last_node())
                            f[node - DoF * first_node()] += integrate_boundary_gradient(be, b, i, quad_shift(b, el), bounds_cond[b].func(comp));
                    }
                }
            }
}

template<class T, class I>
void finite_element_solver_base<T, I>::PETSc_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                    const Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K) {
    Mat A = nullptr;
    MatCreateMPISBAIJWithArrays(PETSC_COMM_WORLD, 1, K.rows(), K.rows(), PETSC_DETERMINE, PETSC_DETERMINE,
                                K.outerIndexPtr(), K.innerIndexPtr(), K.valuePtr(), &A);

    Vec f_petsc = nullptr;
    VecCreate(PETSC_COMM_WORLD, &f_petsc);
    VecSetType(f_petsc, VECSTANDARD);
    VecSetSizes(f_petsc, K.rows(), K.cols());
    for(I i = first_node(); i < last_node(); ++i)
        VecSetValues(f_petsc, 1, &i, &f[i - first_node()], INSERT_VALUES);
    VecAssemblyBegin(f_petsc);
    VecAssemblyEnd(f_petsc);

    Vec x = nullptr;
    VecDuplicate(f_petsc, &x);
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

    KSP ksp = nullptr;
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, f_petsc, x);

    Vec y = nullptr;
    VecScatter toall = nullptr;
    VecScatterCreateToAll(x, &toall, &y);
    VecScatterBegin(toall, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(toall, x, y, INSERT_VALUES, SCATTER_FORWARD);

    f.resize(K.cols());
    PetscScalar* data = nullptr;
    VecGetArray(y, &data);
    for(I i = 0; i < f.size(); ++i)
        f[i] = data[i];

    VecScatterDestroy(&toall);
    KSPDestroy(&ksp);
    MatDestroy(&A);
    VecDestroy(&f_petsc);
    VecDestroy(&x);
    VecDestroy(&y);
}

}

#endif
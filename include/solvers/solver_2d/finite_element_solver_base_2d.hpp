#ifndef FINITE_ELEMENT_ROUTINE_HPP
#define FINITE_ELEMENT_ROUTINE_HPP

#include <numeric>
#ifdef MPI_USE
//#include <mkl.h>
//#include <mkl_cluster_sparse_solver.h>
//#include <petsc.h>
//#include <petscsystypes.h>
#endif
#include <eigen3/Eigen/Sparse>
//#define EIGEN_USE_MKL_ALL
//#include <eigen3/Eigen/PardisoSupport>
#undef I // for new version GCC, when use I macros
#include "mesh.hpp"
#include "right_partition.hpp"
#include "boundary_condition.hpp"

namespace nonlocal {

template<class T, class I, class Matrix_Index>
class finite_element_solver_base {
    std::shared_ptr<mesh::mesh_proxy<T, I>> _mesh_proxy;

protected:
    enum component : bool { X, Y };
    enum theory : bool { LOCAL, NONLOCAL };
    enum index_stage : bool { SHIFTS, NONZERO };

    template<size_t DoF>
    class indexator {
        const std::vector<bool>& _inner_nodes;
        const size_t _node_shift;
        Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& _K_inner,
                                                             & _K_bound;
        std::array<std::vector<bool>, DoF> _inner, _bound;
        std::array<size_t, DoF> _inner_index = {}, _bound_index = {};

    public:
        explicit indexator(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                           Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                           const std::vector<bool>& inner_nodes, const size_t node_shift)
            : _inner_nodes{inner_nodes}
            , _node_shift{node_shift}
            , _K_inner{K_inner}
            , _K_bound{K_bound} {
            for(size_t i = 0; i < DoF; ++i) {
                _inner[i].resize(_inner_nodes.size());
                _bound[i].resize(_inner_nodes.size());
            }
        }

        void fill(const size_t node) {
            for(size_t i = 0; i < DoF; ++i) {
                _inner_index[i] = _K_inner.outerIndexPtr()[DoF * (node - _node_shift) + i];
                _bound_index[i] = _K_bound.outerIndexPtr()[DoF * (node - _node_shift) + i];
                std::fill(_bound[i].begin(), _bound[i].end(), false);
                std::fill(std::next(_inner[i].begin(), DoF * node), _inner[i].end(), false);
            }
        }

        template<index_stage Stage>
        void stage(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K,
                   std::array<std::vector<bool>, DoF>& inner,
                   std::array<size_t, DoF>& inner_index,
                   const size_t row, const size_t col) {
            const size_t i = row % DoF;
            if (!inner[i][col]) {
                if constexpr (Stage == index_stage::SHIFTS)
                    ++K.outerIndexPtr()[row - DoF * _node_shift + 1];
                if constexpr (Stage == index_stage::NONZERO) {
                    K.valuePtr()[inner_index[i]] = 0;
                    K.innerIndexPtr()[inner_index[i]++] = col;
                }
                inner[i][col] = true;
            }
        }

        template<index_stage Stage>
        void index(const size_t block_row, const size_t block_col) {
            for(size_t comp_row = 0; comp_row < DoF; ++comp_row)
                for(size_t comp_col = 0; comp_col < DoF; ++comp_col) {
                    const size_t row = DoF * block_row + comp_row,
                                 col = DoF * block_col + comp_col;
                    if (_inner_nodes[row] && _inner_nodes[col]) {
                        if (row <= col)
                            stage<Stage>(_K_inner, _inner, _inner_index, row, col);
                    } else if (row != col) {
                        if (!_inner_nodes[col])
                            stage<Stage>(_K_bound, _bound, _bound_index, row, col);
                    } else
                        stage<Stage>(_K_inner, _inner, _inner_index, row, col);
                }
        }
    };

    template<size_t DoF, index_stage Stage, theory Theory>
    void mesh_index(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                    const std::vector<bool>& inner_nodes) const {
        indexator<DoF> ind{K_inner, K_bound, inner_nodes, first_node()};
#pragma omp parallel for default(none) firstprivate(ind) schedule(dynamic)
        for(size_t node = first_node(); node < last_node(); ++node) {
            ind.fill(node);
            for(const I eL : mesh_proxy()->nodes_elements_map(node)) {
                if constexpr (Theory == theory::LOCAL)
                    for(size_t jL = 0; jL < mesh().nodes_count(eL); ++jL)
                        ind.template index<Stage>(node, mesh().node_number(eL, jL));
                if constexpr (Theory == theory::NONLOCAL)
                    for(const I eNL : mesh_proxy()->neighbors(eL))
                        for(size_t jNL = 0; jNL < mesh().nodes_count(eNL); ++jNL)
                            ind.template index<Stage>(node, mesh().node_number(eNL, jNL));
            }
        }
    }

    static void prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K) {
        for(size_t i = 0; i < K.rows(); ++i)
            K.outerIndexPtr()[i+1] += K.outerIndexPtr()[i];
        K.data().resize(K.outerIndexPtr()[K.rows()]);
        for(size_t i = 0; i < K.outerIndexPtr()[K.rows()]; ++i)
            K.innerIndexPtr()[i] = K.cols()-1;
    }

    static void sort_indices(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K) {
#pragma omp parallel for default(none) shared(K) schedule(dynamic)
        for(size_t i = 0; i < K.rows(); ++i)
            std::sort(&K.innerIndexPtr()[K.outerIndexPtr()[i]], &K.innerIndexPtr()[K.outerIndexPtr()[i+1]]);
    }

    static constexpr T MAX_LOCAL_WEIGHT = 0.999;

    explicit finite_element_solver_base(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh);
    virtual ~finite_element_solver_base() noexcept = default;

    const mesh::mesh_2d<T, I>&            mesh                     ()                     const { return _mesh_proxy->mesh(); }
    int                                   rank                     ()                     const { return _mesh_proxy->rank(); }
    int                                   size                     ()                     const { return _mesh_proxy->size(); }
    size_t                                first_node               ()                     const { return _mesh_proxy->first_node(); }
    size_t                                last_node                ()                     const { return _mesh_proxy->last_node(); }

    static T jacobian(const std::array<T, 4>& J) noexcept;
    static T jacobian(const std::array<T, 2>& J) noexcept;

    template<class B, size_t DoF>
    std::vector<bool> calc_inner_nodes(const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const;

    template<size_t DoF>
    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                const std::vector<bool>& inner_nodes, const bool nonlocal_task) const;

    template<theory Type, class Callback>
    void mesh_run(const Callback& callback) const;

    template<size_t DoF, class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                     const std::vector<bool>& inner_nodes,
                     const bool nonlocal_task, const Influence_Function& influence_fun,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc) const;

    // Функция обхода групп граничных элементов.
    // Callback - функтор с сигнатурой void(size_t, size_t, size_t)
    template<class Callback>
    void boundary_nodes_run(const Callback& callback) const;

    // Function - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Function>
    T integrate_function(const size_t e, const size_t i, const Function& func) const;

    template<size_t DoF>
    void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                              const right_partition<T, DoF>& right_part) const;

    // Boundary_Gradient - функтор с сигнатурой T(std::array<T, 2>&)
    template<class Boundary_Gradient>
    T integrate_boundary_gradient(const size_t b, const size_t e, const size_t i,
                                  const Boundary_Gradient& boundary_gradient) const;

    template<class B, size_t DoF>
    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const std::vector<boundary_condition<T, B, DoF>>& bounds_cond,
                                       const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound) const;

    template<class B, size_t DoF>
    void integrate_boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                  const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const;

#ifdef MPI_USE
    //void PETSc_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
    //                  const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K);

//    template<size_t DoF>
//    Eigen::Matrix<T, Eigen::Dynamic, 1>
//    MKL_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
//               const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K);
#endif

public:
    void set_mesh(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy);
    const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy() const;
};

template<class T, class I, class Matrix_Index>
finite_element_solver_base<T, I, Matrix_Index>::finite_element_solver_base(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh) { set_mesh(mesh); }

template<class T, class I, class Matrix_Index>
void finite_element_solver_base<T, I, Matrix_Index>::set_mesh(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh_proxy) {
    if (mesh_proxy == nullptr)
        throw std::invalid_argument{"mesh_proxy can't nullptr"};
    _mesh_proxy = mesh_proxy;
}

template<class T, class I, class Matrix_Index>
const std::shared_ptr<mesh::mesh_proxy<T, I>>& finite_element_solver_base<T, I, Matrix_Index>::mesh_proxy() const { return _mesh_proxy; }

template<class T, class I, class Matrix_Index>
T finite_element_solver_base<T, I, Matrix_Index>::jacobian(const std::array<T, 4>& J) noexcept { return std::abs (J[0] * J[3] - J[1] * J[2]); }

template<class T, class I, class Matrix_Index>
T finite_element_solver_base<T, I, Matrix_Index>::jacobian(const std::array<T, 2>& J) noexcept { return std::sqrt(J[0] * J[0] + J[1] * J[1]); }

template<class T, class I, class Matrix_Index>
template<class B, size_t DoF>
std::vector<bool> finite_element_solver_base<T, I, Matrix_Index>::calc_inner_nodes(const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const {
    std::vector<bool> inner_nodes(DoF * mesh().nodes_count(), true);
    boundary_nodes_run(
        [this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < DoF; ++comp)
                if(boundary_type(bounds_cond[b].type(comp)) == boundary_type::FIRST_KIND)
                    inner_nodes[DoF * mesh().node_number(b, el, i) + comp] = false;
        });
    return std::move(inner_nodes);
}

template<class T, class I, class Matrix_Index>
template<size_t DoF>
void finite_element_solver_base<T, I, Matrix_Index>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                                            Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                                            const std::vector<bool>& inner_nodes, const bool nonlocal_task) const {
    if (nonlocal_task) mesh_index<DoF, index_stage::SHIFTS, theory::NONLOCAL>(K_inner, K_bound, inner_nodes);
    else               mesh_index<DoF, index_stage::SHIFTS, theory::   LOCAL>(K_inner, K_bound, inner_nodes);
    prepare_memory(K_inner);
    prepare_memory(K_bound);
    if (nonlocal_task) mesh_index<DoF, index_stage::NONZERO, theory::NONLOCAL>(K_inner, K_bound, inner_nodes);
    else               mesh_index<DoF, index_stage::NONZERO, theory::   LOCAL>(K_inner, K_bound, inner_nodes);
    sort_indices(K_inner);
    sort_indices(K_bound);
}

template<class T, class I, class Matrix_Index>
template<typename finite_element_solver_base<T, I, Matrix_Index>::theory Theory, class Callback>
void finite_element_solver_base<T, I, Matrix_Index>::mesh_run(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback) schedule(dynamic)
    for(size_t node = first_node(); node < last_node(); ++node) {
        for(const I eL : mesh_proxy()->nodes_elements_map(node)) {
            const size_t iL = mesh_proxy()->global_to_local_numbering(eL).find(node)->second; // Проекционные функции
            if constexpr (Theory == theory::LOCAL)
                for(size_t jL = 0; jL < mesh().nodes_count(eL); ++jL) // Аппроксимационные функции
                    callback(eL, iL, jL);
            if constexpr (Theory == theory::NONLOCAL)
                for(const I eNL : mesh_proxy()->neighbors(eL))
                    for(size_t jNL = 0; jNL < mesh().nodes_count(eNL); ++jNL) // Аппроксимационные функции
                        callback(eL, eNL, iL, jNL);
        }
    }
}

template<class T, class I, class Matrix_Index>
template<size_t DoF, class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_solver_base<T, I, Matrix_Index>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                                 Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                                 const std::vector<bool>& inner_nodes,
                                                                 const bool nonlocal_task, const Influence_Function& influence_fun,
                                                                 const Integrate_Loc& integrate_rule_loc,
                                                                 const Integrate_Nonloc& integrate_rule_nonloc) const {
    const auto calc_predicate = [&inner_nodes](const size_t glob_row, const size_t glob_col, const theory flag) {
        bool calc_flag = false;
        for(size_t row = glob_row; !calc_flag && row < glob_row + DoF; ++row)
            for(size_t col = glob_col; !calc_flag && col < glob_col + DoF; ++col)
                calc_flag = inner_nodes[row] && inner_nodes[col] ? row <= col :
                            row != col                           ? !inner_nodes[col] : flag == theory::LOCAL;
        return calc_flag;
    };

    using block_t = std::array<T, DoF * DoF>;
    const auto calc =
        [&K_inner, &K_bound, &inner_nodes, shift = DoF * first_node()]
        (const block_t& block, const size_t glob_row, const size_t glob_col, const theory flag) {
            for(auto [row, component] = std::make_tuple(glob_row, block.cbegin()); row < glob_row + DoF; ++row) {
                T* data_ptr = nullptr;
                for(size_t col = glob_col; col < glob_col + DoF; ++col, ++component)
                    if (inner_nodes[row] && inner_nodes[col]) {
                        if (row <= col) {
                            if (!data_ptr) data_ptr = &K_inner.coeffRef(row - shift, col);
                            else         ++data_ptr;
                            *data_ptr += *component;
                        }
                    } else if (row != col) {
                        if (!inner_nodes[col])
                            K_bound.coeffRef(row - shift, col) += *component;
                    } else if (flag == theory::LOCAL)
                        K_inner.coeffRef(row - shift, col) = T{1};
            }
        };

    mesh_run<theory::LOCAL>(
        [this, &calc_predicate, &calc, &integrate_rule_loc](const size_t e, const size_t i, const size_t j) {
            const I glob_row = DoF * mesh().node_number(e, i),
                    glob_col = DoF * mesh().node_number(e, j);
            if (calc_predicate(glob_row, glob_col, theory::LOCAL))
                calc(block_t{integrate_rule_loc(e, i, j)}, glob_row, glob_col, theory::LOCAL);
        });

    if (nonlocal_task) {
        mesh_run<theory::NONLOCAL>(
            [this, &calc_predicate, &calc, &integrate_rule_nonloc, &influence_fun](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const I glob_row = DoF * mesh().node_number(eL,   iL),
                        glob_col = DoF * mesh().node_number(eNL, jNL);
                if (calc_predicate(glob_row, glob_col, theory::NONLOCAL))
                    calc(block_t{integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun)}, glob_row, glob_col, theory::NONLOCAL);
            });
    }
}

template<class T, class I, class Matrix_Index>
template<class Callback>
void finite_element_solver_base<T, I, Matrix_Index>::boundary_nodes_run(const Callback& callback) const {
    for(size_t b = 0; b < mesh().boundary_groups_count(); ++b)
        for(size_t el = 0; el < mesh().elements_count(b); ++el) {
            const auto& be = mesh().element_1d(b, el);
            for(size_t i = 0; i < be->nodes_count(); ++i)
                callback(b, el, i);
        }
}

template<class T, class I, class Matrix_Index>
template<class Function>
T finite_element_solver_base<T, I, Matrix_Index>::integrate_function(const size_t e, const size_t i, const Function& func) const {
    T integral = 0;
    const auto& el = mesh().element_2d(e);
    auto J = mesh_proxy()->jacobi_matrix(e);
    auto qcoord = mesh_proxy()->quad_coord(e);
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * func(*qcoord) * jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<size_t DoF>
void finite_element_solver_base<T, I, Matrix_Index>::integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                               const right_partition<T, DoF>& right_part) const {
#pragma omp parallel for default(none) shared(f, right_part)
    for(size_t node = first_node(); node < last_node(); ++node)
        for(const I e : mesh_proxy()->nodes_elements_map(node)) {
            const size_t i = mesh_proxy()->global_to_local_numbering(e).find(node)->second;
            for(size_t comp = 0; comp < DoF; ++comp)
                f[DoF * (node - first_node()) + comp] += integrate_function(e, i, right_part[comp]);
        }
}

template<class T, class I, class Matrix_Index>
template<class Boundary_Gradient>
T finite_element_solver_base<T, I, Matrix_Index>::integrate_boundary_gradient(const size_t b, const size_t e, const size_t i,
                                                                              const Boundary_Gradient& boundary_gradient) const {
    T integral = 0;
    const auto& be = mesh().element_1d(b, e);
    auto qcoord = mesh_proxy()->quad_coord(b, e);
    auto J      = mesh_proxy()->jacobi_matrix(b, e);
    for(size_t q = 0; q < be->qnodes_count(); ++q, ++qcoord, ++J)
        integral += be->weight(q) * be->qN(i, q) * boundary_gradient(*qcoord) * jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class B, size_t DoF>
void finite_element_solver_base<T, I, Matrix_Index>::boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                                   const std::vector<boundary_condition<T, B, DoF>>& bounds_cond,
                                                                                   const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound) const {
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

template<class T, class I, class Matrix_Index>
template<class B, size_t DoF>
void finite_element_solver_base<T, I, Matrix_Index>::integrate_boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                                              const std::vector<boundary_condition<T, B, DoF>>& bounds_cond) const {
    for(size_t b = 0; b < bounds_cond.size(); ++b)
        for(size_t comp = 0; comp < DoF; ++comp)
            if(bounds_cond[b].type(comp) == B(boundary_type::SECOND_KIND)) {
                for(size_t e = 0; e < mesh().elements_count(b); ++e) {
                    const auto& be = mesh().element_1d(b, e);
                    for(size_t i = 0; i < be->nodes_count(); ++i) {
                        const I node = DoF * mesh().node_number(b, e, i) + comp;
                        if(node >= DoF * first_node() && node < DoF * last_node())
                            f[node - DoF * first_node()] += integrate_boundary_gradient(b, e, i, bounds_cond[b].func(comp));
                    }
                }
            }
}

#ifdef MPI_USE

/*
template<class T, class I, class Matrix_Index>
void finite_element_solver_base<T, I, Matrix_Index>::PETSc_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                  const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K) {
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
    KSPSetType(ksp, KSPSYMMLQ);
    KSPSetOperators(ksp, A, A);
    KSPSolve(ksp, f_petsc, x);
    int its = 0;
    KSPGetIterationNumber(ksp, &its);
    std::cout << "Iterations = " << its << std::endl;

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
 */

//template<class T, class I, class Matrix_Index>
//template<size_t DoF>
//Eigen::Matrix<T, Eigen::Dynamic, 1>
//finite_element_solver_base<T, I, Matrix_Index>::MKL_solver(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
//                                                           const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K) {
//    void *pt[64] = {};
//    const Matrix_Index maxfct =  1,       // Maximum number of numerical factorizations.
//                       mnum   =  1,       // Which factorization to use.
//                       mtype  =  2,       // real and symmetric positive definite
//                       n      = K.cols(), // matrix size
//                       nrhs   = 1,        // right hands count
//                       msglvl = 1;        // Print statistical information in file
//    T fdum = 0;                           // ignored
//    Matrix_Index idum = 0;                // ignored
//    Matrix_Index error = 0;               // error flag
//    MPI_Fint comm = MPI_Comm_c2f(MPI_COMM_WORLD);
//
//    // https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/parallel-direct-sparse-solver-for-clusters-interface/cluster-sparse-solver-iparm-parameter.html#cluster-sparse-solver-iparm-parameter
//    std::array<Matrix_Index, 64> iparm = {};
//    iparm[ 0] =  1; // Solver default parameters overriden with provided by iparm
//    iparm[ 1] =  3; // reordering
//    iparm[ 7] =  5; // Max number of iterative refinement steps
//    iparm[ 9] = 13; // Perturb the pivot elements with 1E-13
//    iparm[17] = -1; // Output: Number of nonzeros in the factor LU
//    iparm[18] = -1; // Output: Mflops for LU factorization
//    iparm[27] = std::is_same_v<T, float>;
//    iparm[34] =  1; // Cluster Sparse Solver use C-style indexing for ia and ja arrays
//    iparm[39] =  3; // Input: matrix/rhs are distributed between MPI processes
//    iparm[40] = DoF * first_node();
//    iparm[41] = iparm[40] + K.rows() - 1;
//
//    Eigen::Matrix<T, Eigen::Dynamic, 1> x = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(K.cols());
//    if constexpr (std::is_same_v<Matrix_Index, int>) {
//        Matrix_Index phase = 11;
//        cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                              &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//        phase = 22;
//        cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                              &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//        phase = 33;
//        cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                              &idum, &nrhs, iparm.data(), &msglvl, f.data(), x.data(), &comm, &error);
//        phase = -1; // Release internal memory.
//        cluster_sparse_solver(pt, &maxfct, &mnum, &mtype, &phase, &n, &fdum, K.outerIndexPtr(), K.innerIndexPtr(),
//                              &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//    } else {
//        Matrix_Index phase = 11;
//        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                                 &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//        phase = 22;
//        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                                 &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//        phase = 33;
//        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &n, K.valuePtr(), K.outerIndexPtr(), K.innerIndexPtr(),
//                                 &idum, &nrhs, iparm.data(), &msglvl, f.data(), x.data(), &comm, &error);
//        phase = -1; // Release internal memory.
//        cluster_sparse_solver_64(pt, &maxfct, &mnum, &mtype, &phase, &n, &fdum, K.outerIndexPtr(), K.innerIndexPtr(),
//                                 &idum, &nrhs, iparm.data(), &msglvl, &fdum, &fdum, &comm, &error);
//    }
//
//    static constexpr std::array<std::string_view, 12> errors = {
//        "no error",
//        "input inconsistent",
//        "not enough memory",
//        "reordering problem",
//        "Zero pivot, numerical factorization or iterative refinement problem. If the error appears during the solution phase, try to change the pivoting perturbation (iparm[9]) and also increase the number of iterative refinement steps. If it does not help, consider changing the scaling, matching and pivoting options (iparm[10], iparm[12], iparm[20])",
//        "unclassified (internal) error",
//        "reordering failed (matrix types 11 and 13 only)",
//        "diagonal matrix is singular",
//        "32-bit integer overflow problem",
//        "not enough memory for OOC",
//        "error opening OOC files",
//        "read/write error with OOC files"
//    };
//    error = std::abs(error);
//    if (error < 13)
//        std::cerr << errors[error] << std::endl;
//
//    return std::move(x);
//}
#endif

}

#endif
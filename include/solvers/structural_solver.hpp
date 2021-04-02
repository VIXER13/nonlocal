#ifndef STRUCTURAL_SOLVER_HPP
#define STRUCTURAL_SOLVER_HPP

#include <functional>
#include <algorithm>
#include <omp.h>
#include "finite_element_solver_base.hpp"
#include "structural_solution.hpp"
#include "../../Eigen/Eigen/Dense"
#include "../../Eigen/Eigen/Sparse"
//#include "Eigen/PardisoSupport"

namespace nonlocal::structural {

enum class boundary_t : uint8_t {
    DISPLACEMENT = uint8_t(boundary_type::FIRST_KIND),
    PRESSURE     = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
using bound_cond = boundary_condition<T, boundary_t, 2>;

template<class T>
using right_part = right_partition<T, 2>;

template<class T, class I>
class structural_solver : public finite_element_solver_base<T, I> {
    using _base = finite_element_solver_base<T, I>;
    using typename _base::component;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::jacobian;
    using _base::first_node;
    using _base::last_node;

    std::array<T, 3> _D;

    template<bool Proj, bool Approx>
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;

    template<bool Proj, bool Approx, class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                const std::vector<bool>& inner_nodes,
                                const bool nonlocal_task) const;

    template<class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                     const T p1, const Influence_Function& influence_fun,
                     const std::vector<bool>& inner_nodes) const;

    template<class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                       Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const T p1, const Influence_Function& influence_fun);

    template<class Influence_Function>
    void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                               const T alpha, const std::vector<T>& temperature,
                               const T p1, const Influence_Function& influence_fun);

public:
    explicit structural_solver(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh)
        : _base{mesh} {}

    ~structural_solver() override = default;

    template<class Influence_Function>
    solution<T, I> stationary(const calculation_parameters<T>& params,
                              const std::vector<bound_cond<T>> &bounds_cond,
                              const right_part<T>& right_part,
                              const T p1, const Influence_Function& influence_fun);

    //template<class Influence_Function>
    //solution<T, I> stationary(const std::vector<bound_cond<T>> &bounds_cond,
    //                          const right_part<T>& right_part,
    //                          const T alpha, const std::vector<T>& temperature,
    //                          const T p1, const Influence_Function& influence_fun);
};

template<class T, class I>
template<bool Proj, bool Approx>
T structural_solver<T, I>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    static constexpr size_t k = Proj ^ Approx;
    T integral = 0;
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj)
        integral += el->weight(q) / jacobian(*J) *
                    (_D[k] * (*dNdi)[ Proj] * (*dNdj)[ Approx] +
                     _D[2] * (*dNdi)[!Proj] * (*dNdj)[!Approx]);
    return integral;
}

template<class T, class I>
template<bool Proj, bool Approx, class Influence_Function>
T structural_solver<T, I>::integrate_nonloc(const size_t eL, const size_t eNL,
                                            const size_t iL, const size_t jNL,
                                            const Influence_Function& influence_function) const {
    static constexpr size_t k = Proj ^ Approx;
    T integral = 0;
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_start = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_start    = _base::mesh_proxy()->dNdX(eNL, jNL);
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        std::array<T, 2> inner_int = {};
        auto qcoordNL = qcoordNL_start;
        auto dNdNL    = dNdNL_start;
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_int[X] += influence_weight * (*dNdNL)[ Approx];
            inner_int[Y] += influence_weight * (*dNdNL)[!Approx];
        }
        integral += elL->weight(qL) * (_D[k] * inner_int[X] * (*dNdL)[ Proj] +
                                       _D[2] * inner_int[Y] * (*dNdL)[!Proj]);
    }
    return integral;
}

template<class T, class I>
void structural_solver<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                                                     Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                                     const std::vector<bool>& inner_nodes,
                                                     const bool nonlocal_task) const {
    if (nonlocal_task)
        _base::template mesh_index<2, _base::index_stage::SHIFTS, _base::theory::NONLOCAL>(K_inner, K_bound, inner_nodes);
    else
        _base::template mesh_index<2, _base::index_stage::SHIFTS, _base::theory::LOCAL>(K_inner, K_bound, inner_nodes);

    _base::prepare_memory(K_inner);
    _base::prepare_memory(K_bound);

    if (nonlocal_task)
        _base::template mesh_index<2, _base::index_stage::NONZERO, _base::theory::NONLOCAL>(K_inner, K_bound, inner_nodes);
    else
        _base::template mesh_index<2, _base::index_stage::NONZERO, _base::theory::LOCAL>(K_inner, K_bound, inner_nodes);

    _base::sort_indices(K_inner);
    _base::sort_indices(K_bound);
}

template<class T, class I>
template<class Influence_Function>
void structural_solver<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                          Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                          const T p1, const Influence_Function& influence_fun,
                                          const std::vector<bool>& inner_nodes) const {
    const auto filler_loc =
        [this, &K, &K_bound, &inner_nodes, p1, shift = 2 * first_node()]
        (const size_t e, const size_t i, const size_t j, const component proj, const component approx, const auto& integrate_rule) {
            const I row = 2 * mesh().node_number(e, i) + proj,
                    col = 2 * mesh().node_number(e, j) + approx;
            if (inner_nodes[row] && inner_nodes[col]) {
                if (row <= col)
                    K.coeffRef(row - shift, col) += p1 * integrate_rule(e, i, j);
            } else if (row != col) {
                if (!inner_nodes[col])
                    K_bound.coeffRef(row - shift, col) += p1 * integrate_rule(e, i, j);
            } else
                K.coeffRef(row - shift, col) = 1;
        };

    _base::template mesh_run<_base::theory::LOCAL>(
        [this, &filler_loc](const size_t e, const size_t i, const size_t j) {
            filler_loc(e, i, j, X, X, [this](const size_t e, const size_t i, const size_t j) { return integrate_loc<X, X>(e, i, j); });
            filler_loc(e, i, j, X, Y, [this](const size_t e, const size_t i, const size_t j) { return integrate_loc<X, Y>(e, i, j); });
            filler_loc(e, i, j, Y, X, [this](const size_t e, const size_t i, const size_t j) { return integrate_loc<Y, X>(e, i, j); });
            filler_loc(e, i, j, Y, Y, [this](const size_t e, const size_t i, const size_t j) { return integrate_loc<Y, Y>(e, i, j); });
        });

    if (p1 < _base::MAX_LOCAL_WEIGHT) {
        const auto filler_nonloc =
            [this, &K, &K_bound, &inner_nodes, &influence_fun, p2 = 1 - p1, shift = 2 * first_node()]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const component proj, const component approx, const auto& integrate_rule) {
                const I row = 2 * mesh().node_number(eL,  iL ) + proj,
                        col = 2 * mesh().node_number(eNL, jNL) + approx;
                if (inner_nodes[row] && inner_nodes[col]) {
                    if (row <= col)
                        K.coeffRef(row - shift, col) += p2 * integrate_rule(eL, eNL, iL, jNL, influence_fun);
                } else if (row != col)
                    if (!inner_nodes[col])
                        K_bound.coeffRef(row - shift, col) += p2 * integrate_rule(eL, eNL, iL, jNL, influence_fun);
            };

        _base::template mesh_run<_base::theory::NONLOCAL>(
            [this, &filler_nonloc](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
#define INTEGRATE_NONLOC(X, Y) filler_nonloc(eL, eNL, iL, jNL, X, Y, \
                               [this](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) \
                               { return integrate_nonloc<X, Y, Influence_Function>(eL, eNL, iL, jNL, influence_function); });
                INTEGRATE_NONLOC(X, X)
                INTEGRATE_NONLOC(X, Y)
                INTEGRATE_NONLOC(Y, X)
                INTEGRATE_NONLOC(Y, Y)
#undef INTEGRATE_NONLOC
            });
    }
}

template<class T, class I>
template<class Influence_Function>
void structural_solver<T, I>::create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                            Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                            const std::vector<bound_cond<T>>& bounds_cond,
                                            const T p1, const Influence_Function& influence_fun) {
    std::vector<bool> inner_nodes(2*mesh().nodes_count(), true);
    _base::template boundary_nodes_run(
        [this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
            for(size_t comp = 0; comp < 2; ++comp)
                if(bounds_cond[b].type(comp) == boundary_t::DISPLACEMENT)
                    inner_nodes[2 * mesh().node_number(b, el, i) + comp] = false;
        });

    double time = omp_get_wtime();
    create_matrix_portrait(K, K_bound, inner_nodes, p1 < _base::MAX_LOCAL_WEIGHT);
    std::cout << "create_matrix_portrait: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    calc_matrix(K, K_bound, p1, influence_fun, inner_nodes);
    std::cout << "calc coeffs: " << omp_get_wtime() - time << std::endl;
}

template<class T, class I>
template<class Influence_Function>
solution<T, I> structural_solver<T, I>::stationary(const calculation_parameters<T>& parameters,
                                                   const std::vector<bound_cond<T>> &bounds_cond,
                                                   const right_part<T>& right_part,
                                                   const T p1, const Influence_Function& influence_fun) {
    _D = hooke_matrix(parameters.nu, parameters.E, parameters.type);
    double time = omp_get_wtime();
    const size_t rows = 2 * (last_node() - first_node()),
                 cols = 2 * mesh().nodes_count();
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K      (rows, cols),
                                               K_bound(rows, cols);
    create_matrix(K, K_bound, bounds_cond, p1, influence_fun);
    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    _base::template integrate_right_part(f, right_part);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
    //Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
    _base::PETSc_solver(f, K);
    Eigen::Matrix<T, Eigen::Dynamic, 1> u = f;
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

    return solution<T, I>{_base::mesh_proxy(), parameters, p1, influence_fun, u};
}

//template<class T, class I>
//template<class Influence_Function>
//void structural_solver<T, I>::temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
//                                                    const T alpha, const std::vector<T>& temperature,
//                                                    const T p1, const Influence_Function& influence_fun) {
////    const T c = alpha * _params.E / (1 - 2 * _params.nu);
////    const auto [Tx, Ty] = _base::get_mesh_proxy()->calc_gradient(temperature);
////    const std::vector<T> qTx = _base::get_mesh_proxy()->approx_in_quad(Tx),
////                         qTy = _base::get_mesh_proxy()->approx_in_quad(Ty);
////
////#pragma omp parallel for default(none) shared(f)
////    for(size_t node = 0; node < mesh().nodes_count(); ++node) {
////        for(const I e : _base::nodes_elements_map(node)) {
////            const auto& el = mesh().element_2d(e);
////            const size_t i = _base::global_to_local_numbering(e).find(node)->second;
////            T int_x = 0, int_y = 0;
////            for(size_t q = 0, shift = quad_shift(e); q < el->qnodes_count(); ++q, ++shift) {
////                const T jac = el->weight(q) * jacobian(shift) * el->qN(i, q);
////                int_x += qTx[shift] * jac;
////                int_y += qTy[shift] * jac;
////            }
////            f[2 * node]   += p1 * c * int_x;
////            f[2 * node+1] += p1 * c * int_y;
////        }
////    }
////
////    if(p1 < _base::MAX_LOCAL_WEIGHT) {
////        const T p2 = 1 - p1;
////#pragma omp parallel for default(none) shared(f) firstprivate(influence_fun)
////        for(size_t node = 0; node < mesh().nodes_count(); ++node) {
////            for(const I eL : _base::nodes_elements_map(node)) {
////                const auto& elL = mesh().element_2d(eL);
////                const size_t iL = _base::global_to_local_numbering(eL).find(node)->second;
////                T int_x = 0, int_y = 0;
////                for(size_t qL = 0, shiftL = quad_shift(eL); qL < elL->qnodes_count(); ++qL, ++shiftL) {
////                    const T jac = elL->weight(qL) * jacobian(shiftL) * elL->qN(iL, qL);
////                    T inner_int_x = 0, inner_int_y = 0;
////                    for(const I eNL : _base::neighbors(eL)) {
////                        const auto& elNL = mesh().element_2d(eNL);
////                        for(size_t qNL = 0, shiftNL = quad_shift(eNL); qNL < elNL->qnodes_count(); ++qNL, ++shiftNL) {
////                            const T influence_weight = elNL->weight(qNL) * jacobian(shiftNL) * influence_fun(quad_coord(shiftL), quad_coord(shiftNL));
////                            inner_int_x += influence_weight * qTx[shiftNL];
////                            inner_int_y += influence_weight * qTy[shiftNL];
////                        }
////                    }
////                    int_x += inner_int_x * jac;
////                    int_y += inner_int_y * jac;
////                }
////                f[2 * node]   += p1 * c * int_x;
////                f[2 * node+1] += p1 * c * int_y;
////            }
////        }
////    }
//}
//
//template<class T, class I>
//template<class Influence_Function>
//solution<T, I> structural_solver<T, I>::stationary(const std::vector<bound_cond<T>> &bounds_cond,
//                                                   const right_part<T>& right_part,
//                                                   const T alpha, const std::vector<T>& temperature,
//                                                   const T p1, const Influence_Function& influence_fun) {
//    double time = omp_get_wtime();
//    const size_t rows = 2 * (last_node() - first_node()),
//                 cols = 2 * mesh().nodes_count();
//    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K      (rows, cols),
//                                               K_bound(rows, cols);
//    create_matrix(K, K_bound, bounds_cond, p1, influence_fun);
//    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;
//
//
//    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
//
//    time = omp_get_wtime();
//    _base::template integrate_right_part(f, right_part);
//    std::cout << "f: " << omp_get_wtime() - time << std::endl;
//
//    time = omp_get_wtime();
//    temperature_condition(f, alpha, temperature, p1, influence_fun);
//    std::cout << "temp cond: " << omp_get_wtime() - time << std::endl;
//
//    time = omp_get_wtime();
//    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
//    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
//    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;
//
//    time = omp_get_wtime();
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
//    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
//    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;
//
//    return solution<T, I>{_base::get_mesh_proxy(), _D, u, p1, influence_fun};
//}

}

#endif
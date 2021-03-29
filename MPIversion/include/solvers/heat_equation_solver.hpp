#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <iostream>
#include <algorithm>
#include <omp.h>
#include "finite_element_solver_base.hpp"
#include "heat_equation_solution.hpp"
#include "../../Eigen/Eigen/Dense"
#include "../../Eigen/Eigen/Sparse"
#include "../../Eigen/Eigen/Eigen"
#include "../../Eigen/Eigen/SparseCholesky"
//#include "Eigen/PardisoSupport"

namespace nonlocal::heat {

enum class boundary_t : uint8_t {
    TEMPERATURE = uint8_t(boundary_type::FIRST_KIND),
    FLOW        = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
using bound_cond = boundary_condition<T, boundary_t, 1>;

template<class T>
using right_part = right_partition<T, 1>;

template<class T, class I>
class heat_equation_solver : public finite_element_solver_base<T, I> {
    using _base = finite_element_solver_base<T, I>;
    using _base::size;
    using _base::rank;
    using _base::first_node;
    using _base::last_node;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::jacobian;

    static constexpr size_t DoF = 1;

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                const bool neumann_task, const bool nonlocal_task,
                                const std::vector<bool>& inner_nodes) const;

    template<class Integrate_Rule, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                     const Integrate_Rule& integrate_rule, const bool neumann_task,
                     const T p1, const Influence_Function& influence_fun,
                     const std::vector<bool>& inner_nodes) const;

    template<class Integrate_Rule, class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                       Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const bool neumann_task, const Integrate_Rule& integrate_rule,
                       const T p1, const Influence_Function& influence_fun) const;

public:
    explicit heat_equation_solver(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh)
        : _base{mesh} {}

    ~heat_equation_solver() override = default;

    template<class Influence_Function>
    solution<T, I> stationary(const std::vector<bound_cond<T>>& bounds_cond, const right_part<T>& right_part,
                              const T p1, const Influence_Function& influence_fun, const T volume = 0);

//    template<class Init_Distribution, class Right_Part, class Influence_Function>
//    void nonstationary(const std::string& path, const T tau, const uintmax_t time_steps,
//                       const std::vector<bound_cond<T>>& bounds_cond,
//                       const Init_Distribution& init_dist, const Right_Part& right_part,
//                       const T r, const T p1, const Influence_Function& influence_fun,
//                       const uintmax_t print_frequency = std::numeric_limits<uintmax_t>::max());
};

template<class T, class I>
T heat_equation_solver<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * jacobian(*J);
    return integral;
}

template<class T, class I>
T heat_equation_solver<T, I>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q) * jacobian(*J);
    return integral;
}

template<class T, class I>
T heat_equation_solver<T, I>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj)
        integral += el->weight(q) * ((*dNdi)[X] * (*dNdj)[X] + (*dNdi)[Y] * (*dNdj)[Y]) / jacobian(*J);
    return integral;
}

template<class T, class I>
template<class Influence_Function>
T heat_equation_solver<T, I>::integrate_nonloc(const size_t eL, const size_t eNL,
                                               const size_t iL, const size_t jNL,
                                               const Influence_Function& influence_function) const {
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
            inner_int[X] += influence_weight * (*dNdNL)[X];
            inner_int[Y] += influence_weight * (*dNdNL)[Y];
        }
        integral += elL->weight(qL) * (inner_int[X] * (*dNdL)[X] + inner_int[Y] * (*dNdL)[Y]);
    }
    return integral;
}

template<class T, class I>
void heat_equation_solver<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_inner,
                                                        Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                                        const bool neumann_task, const bool nonlocal_task,
                                                        const std::vector<bool>& inner_nodes) const {
    const auto outer_index_filler =
        [&K_inner, &K_bound, neumann_task, shift = I(first_node())-1]
        (const size_t node, const std::vector<bool>& inner, const std::vector<bool>& bound) {
            K_inner.outerIndexPtr()[node - shift] = std::accumulate(std::next(inner.cbegin(), node), inner.cend(), size_t{0}) + neumann_task;
            K_bound.outerIndexPtr()[node - shift] = std::accumulate(bound.cbegin(), bound.cend(), size_t{0});
        };

    if (nonlocal_task)
        _base::template mesh_index<DoF, _base::task_type::NONLOCAL>(inner_nodes, outer_index_filler);
    else
        _base::template mesh_index<DoF, _base::task_type::LOCAL>(inner_nodes, outer_index_filler);

    for(size_t i = 0; i < K_inner.rows(); ++i)
        K_inner.outerIndexPtr()[i+1] += K_inner.outerIndexPtr()[i];
    K_inner.data().resize(K_inner.outerIndexPtr()[K_inner.rows()]);

    for(size_t i = 0; i < K_bound.rows(); ++i)
        K_bound.outerIndexPtr()[i+1] += K_bound.outerIndexPtr()[i];
    K_bound.data().resize(K_bound.outerIndexPtr()[K_bound.rows()]);

    const auto inner_index_filler =
        [&K_inner, &K_bound, neumann_task, shift = first_node()]
        (const size_t node, const std::vector<bool>& inner, const std::vector<bool>& bound) {
            size_t inner_curr = K_inner.outerIndexPtr()[node - shift];
            for(size_t k = node; k < inner.size(); ++k)
                if (inner[k]) {
                    K_inner.valuePtr()[inner_curr] = 0;
                    K_inner.innerIndexPtr()[inner_curr++] = k;
                }

            if (neumann_task) {
                K_inner.valuePtr()[inner_curr] = 0;
                K_inner.innerIndexPtr()[inner_curr] = inner.size();
            }

            for(size_t k = 0, bound_curr = K_bound.outerIndexPtr()[node - shift]; k < bound.size(); ++k)
                if (bound[k]) {
                    K_bound.valuePtr()[bound_curr] = 0;
                    K_bound.innerIndexPtr()[bound_curr++] = k;
                }
        };

    if (nonlocal_task)
        _base::template mesh_index<DoF, _base::task_type::NONLOCAL>(inner_nodes, inner_index_filler);
    else
        _base::template mesh_index<DoF, _base::task_type::LOCAL>(inner_nodes, inner_index_filler);
}

//    template<class T, class I>
//    void heat_equation_solver<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
//                                                            Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
//                                                            const bool neumann_task, const T p1,
//                                                            const std::vector<bool>& inner_nodes) const {
//        std::vector<std::unordered_set<I>> inner_portrait(K.rows()),
//                                           bound_portrait(K_bound.rows());
//        if (neumann_task) {
//#pragma omp parallel for default(none) shared(K, inner_portrait)
//            for(size_t node =  0; node < inner_portrait.size(); ++node)
//                inner_portrait[node].insert(mesh().nodes_count());
//        }
//
//        const auto indexator = [&inner_nodes, &inner_portrait, &bound_portrait, shift = first_node()](const I row, const I col) {
//            if (inner_nodes[row] && inner_nodes[col]) {
//                if (row <= col)
//                    inner_portrait[row - shift].insert(col);
//            } else if (row != col) {
//                if (!inner_nodes[col])
//                    bound_portrait[row - shift].insert(col);
//            } else
//                inner_portrait[row - shift].insert(col);
//        };
//
//        if (p1 > _base::MAX_LOCAL_WEIGHT) {
//            _base::template mesh_run_loc(
//                    [this, &indexator] (const size_t e, const size_t i, const size_t j) {
//                        indexator(mesh().node_number(e, i), mesh().node_number(e, j));
//                    });
//        } else {
//            _base::template mesh_run_nonloc(
//                    [this, &indexator](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
//                        indexator(mesh().node_number(eL, iL), mesh().node_number(eNL, jNL));
//                    });
//        }
//
//        _base::convert_portrait(K,       inner_portrait);
//        _base::convert_portrait(K_bound, bound_portrait);
//    }

template<class T, class I>
template<class Integrate_Rule, class Influence_Function>
void heat_equation_solver<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                             Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                             const Integrate_Rule& integrate_rule, const bool neumann_task,
                                             const T p1, const Influence_Function& influence_fun,
                                             const std::vector<bool>& inner_nodes) const {
    if (neumann_task) {
#pragma omp parallel for default(none) shared(K)
        for(size_t node = first_node(); node < last_node(); ++node) {
            T& val = K.coeffRef(node - first_node(), mesh().nodes_count());
            for(const I e : _base::mesh_proxy()->nodes_elements_map(node))
                val += integrate_basic(e, _base::mesh_proxy()->global_to_local_numbering(e).find(node)->second);
        }
    }

    _base::template mesh_run_loc(
        [this, &K, &K_bound, &inner_nodes, &integrate_rule, p1](const size_t e, const size_t i, const size_t j) {
            const I row = mesh().node_number(e, i),
                    col = mesh().node_number(e, j);
            if (inner_nodes[row] && inner_nodes[col]) {
                if (row <= col)
                    K.coeffRef(row - first_node(), col) += p1 * integrate_rule(e, i, j);
            } else if (row != col) {
                if (!inner_nodes[col])
                    K_bound.coeffRef(row - first_node(), col) += p1 * integrate_rule(e, i, j);
            } else
                K.coeffRef(row - first_node(), col) = 1;
        }
    );

    if (p1 < _base::MAX_LOCAL_WEIGHT) {
        _base::template mesh_run_nonloc(
            [this, &K, &K_bound, &inner_nodes, &influence_fun, p2 = 1 - p1]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const I row = mesh().node_number(eL,  iL ),
                        col = mesh().node_number(eNL, jNL);
                if (inner_nodes[row] && inner_nodes[col]) {
                    if (row <= col)
                        K.coeffRef(row - first_node(), col) += p2 * integrate_nonloc(eL, eNL, iL, jNL, influence_fun);
                } else if (row != col)
                    if (!inner_nodes[col])
                        K_bound.coeffRef(row - first_node(), col) += p2 * integrate_nonloc(eL, eNL, iL, jNL, influence_fun);
            }
        );
    }
}

template<class T, class I>
template<class Integrate_Rule, class Influence_Function>
void heat_equation_solver<T, I>::create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                               Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                               const std::vector<bound_cond<T>>& bounds_cond,
                                               const bool neumann_task, const Integrate_Rule& integrate_rule,
                                               const T p1, const Influence_Function& influence_fun) const {
    std::vector<bool> inner_nodes(mesh().nodes_count(), true);
    _base::template boundary_nodes_run([this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
        if(bounds_cond[b].type(0) == boundary_t::TEMPERATURE)
            inner_nodes[mesh().node_number(b, el, i)] = false;
    });

    double time = omp_get_wtime();
    create_matrix_portrait(K, K_bound, neumann_task, p1 < _base::MAX_LOCAL_WEIGHT, inner_nodes);
    std::cout << "rank = " << rank() << std::endl;
    std::cout << "K.nonzero() = " << K.nonZeros() << std::endl;
    std::cout << "K_bound.nonzero() = " << K_bound.nonZeros() << std::endl;
    std::cout << "create_matrix_portrait: " << omp_get_wtime() - time << std::endl << std::endl;

    time = omp_get_wtime();
    calc_matrix(K, K_bound, integrate_rule, neumann_task, p1, influence_fun, inner_nodes);
    std::cout << "rank = " << rank() << std::endl;
    std::cout << "calc coeffs: " << omp_get_wtime() - time << std::endl;
}

template<class T, class I>
template<class Influence_Function>
solution<T, I> heat_equation_solver<T, I>::stationary(const std::vector<bound_cond<T>>& bounds_cond, const right_part<T>& right_part,
                                                      const T p1, const Influence_Function& influence_fun, const T volume) {
    const bool neumann_task = std::all_of(bounds_cond.cbegin(), bounds_cond.cend(),
                                          [](const bound_cond<T>& bound) { return bound.type(0) == boundary_t::FLOW; });
    const size_t rows = last_node() - first_node() + (neumann_task && rank() == size() - 1),
                 cols = mesh().nodes_count() + neumann_task;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
//        if(neumann_task) {
//            if(std::abs(std::accumulate(f.data(), f.data() + f.size(), T{0}, [](const T sum, const T val) { return sum + val; })) > 1e-5)
//                throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
//            f[mesh().nodes_count()] = volume;
//        }

    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K      (rows, cols),
                                               K_bound(rows, cols);
    create_matrix(
        K, K_bound, bounds_cond, neumann_task,
        [this](const size_t e, const size_t i, const size_t j) { return integrate_loc(e, i, j); },
        p1, influence_fun
    );

    _base::template integrate_right_part(f, right_part);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);

    double time = omp_get_wtime();
    _base::PETSc_solver(f, K);
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

//    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
//
//    //Eigen::PardisoLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::Lower> solver;
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
//    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
//    //std::cout << "System solving: " << omp_get_wtime() - time << std::endl;
//    temperature.conservativeResize(mesh().nodes_count());
//    return std::move(temperature);

    return solution<T, I>{_base::mesh_proxy(), f};
}

//template<class T, class I>
//template<class Init_Distribution, class Right_Part, class Influence_Function>
//void heat_equation_solver<T, I>::nonstationary(const std::string& path, const T tau, const uintmax_t time_steps,
//                                               const std::vector<bound_cond<T>>& bounds_cond,
//                                               const Init_Distribution& init_dist, const Right_Part& right_part,
//                                               const T r, const T p1, const Influence_Function& influence_fun, const uintmax_t print_frequency) {
//    static constexpr bool NOT_NEUMANN_TASK = false;
//    Eigen::SparseMatrix<T, Eigen::ColMajor, I> K      (mesh().nodes_count(), mesh().nodes_count()),
//                                               K_bound(mesh().nodes_count(), mesh().nodes_count());
//    create_matrix(
//       K, K_bound, bounds_cond, NOT_NEUMANN_TASK,
//       [this](const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) {
//            return integrate_loc(e, i, j, quad_shift); },
//       p1, influence_fun
//    );
//
//    static constexpr T LOCAL = 1;
//    Eigen::SparseMatrix<T, Eigen::ColMajor, I> C      (mesh().nodes_count(), mesh().nodes_count()),
//                                               C_bound(mesh().nodes_count(), mesh().nodes_count());
//    create_matrix(
//       C, C_bound, bounds_cond, NOT_NEUMANN_TASK,
//       [this](const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) {
//            return integrate_basic_pair(e, i, j, quad_shift); },
//       LOCAL, influence_fun
//    );
//
//    C_bound.setZero();
//    K_bound *= tau;
//    K *= tau;
//    K += C;
//    _base::template boundary_nodes_run([this, &bounds_cond, &K](const size_t b, const size_t el, const size_t i) {
//        if(bounds_cond[b].type == boundary_t::TEMPERATURE)
//            K.coeffRef(mesh().node_number(b, el, i), mesh().node_number(b, el, i)) = 1;
//    });
//
//    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh().nodes_count()),
//                                        temperature_prev(mesh().nodes_count()),
//                                        temperature     (mesh().nodes_count());
//    for(size_t i = 0; i < mesh().nodes_count(); ++i)
//        temperature_prev[i] = init_dist(mesh().node(i));
//
////    Eigen::PardisoLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::Lower> solver;
////    solver.compute(K);
////    if(print_frequency != std::numeric_limits<uintmax_t>::max()) {
////        save_as_vtk(path + "0.vtk", temperature_prev);
////        std::cout << "step = " << 0 << " Volume = " << integrate_solution(temperature_prev) << std::endl;
////    }
////    for(size_t i = 1; i < time_steps; ++i) {
////        f.setZero();
////        integrate_boundary_flow(f, bounds_cond);
////        integrate_right_part(f, right_part);
////        f *= tau;
////        f += C.template selfadjointView<Eigen::Lower>() * temperature_prev;
////        temperature_on_boundary(f, bounds_cond, K_bound);
////        temperature = solver.solve(f);
////        temperature_prev.swap(temperature);
////        if(i % print_frequency == 0) {
////            save_as_vtk(path + std::to_string(i) + ".vtk", temperature_prev);
////            std::cout << "step = " << i << " Volume = " << integrate_solution(temperature_prev) << std::endl;
////        }
////    }
//}

}

#endif
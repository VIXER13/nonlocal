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
    using typename _base::Finite_Element_2D_Ptr;
    using typename _base::component;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::quad_shift;
    using _base::quad_coord;
    using _base::jacobian;
    using _base::first_node;
    using _base::last_node;

    parameters<T> _params;
    std::array<T, 3> _D;

    template<bool Proj, bool Approx>
    T integrate_loc(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) const;

    template<bool Proj, bool Approx, class Influence_Function>
    T integrate_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,  size_t shiftL,
                       const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                const T p1, const std::vector<bool>& inner_nodes) const;

    template<class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                     const T p1, const Influence_Function& influence_fun,
                     const std::vector<bool>& inner_nodes) const;

    template<class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                       Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const T p1, const Influence_Function& influence_fun);

    template<class Influence_Function>
    void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                               const T alpha, const std::vector<T>& temperature,
                               const T p1, const Influence_Function& influence_fun);

public:
    explicit structural_solver(const std::shared_ptr<mesh::mesh_info<T, I>>& mesh, const parameters<T>& params)
        : _base{mesh}
        , _params{params}
        , _D{hooke_matrix(params)} {}

    ~structural_solver() override = default;

    template<class Influence_Function>
    solution<T, I> stationary(const std::vector<bound_cond<T>> &bounds_cond,
                              const right_part<T>& right_part,
                              const T p1, const Influence_Function& influence_fun);

    template<class Influence_Function>
    solution<T, I> stationary(const std::vector<bound_cond<T>> &bounds_cond,
                              const right_part<T>& right_part,
                              const T alpha, const std::vector<T>& temperature,
                              const T p1, const Influence_Function& influence_fun);
};

template<class T, class I>
template<bool Proj, bool Approx>
T structural_solver<T, I>::integrate_loc(const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift) const {
    T integral = 0;
    static constexpr size_t k = Proj ^ Approx;
    for(size_t q = 0; q < e->qnodes_count(); ++q, ++quad_shift)
        integral += e->weight(q) / jacobian(quad_shift) *
                    (_D[k] * _base::template dNd< Proj>(e, i, q, quad_shift) * _base::template dNd< Approx>(e, j, q, quad_shift) +
                     _D[2] * _base::template dNd<!Proj>(e, i, q, quad_shift) * _base::template dNd<!Approx>(e, j, q, quad_shift));
    return integral;
}

template<class T, class I>
template<bool Proj, bool Approx, class Influence_Function>
T structural_solver<T, I>::integrate_nonloc(const Finite_Element_2D_Ptr& eL,  const size_t iL,  size_t shiftL,
                                            const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL,
                                            const Influence_Function& influence_function) const {
    T integral = 0;
    const size_t sub_shift = shiftNL;
    static constexpr size_t k = Proj ^ Approx;
    for(size_t qL = 0; qL < eL->qnodes_count(); ++qL, ++shiftL) {
        T inner_int_x = 0, inner_int_y = 0;
        for(size_t qNL = 0, shiftNL = sub_shift; qNL < eNL->qnodes_count(); ++qNL, ++shiftNL) {
            const T influence_weight = eNL->weight(qNL) * influence_function(quad_coord(shiftL), quad_coord(shiftNL));
            inner_int_x += influence_weight * _base::template dNd< Approx>(eNL, jNL, qNL, shiftNL);
            inner_int_y += influence_weight * _base::template dNd<!Approx>(eNL, jNL, qNL, shiftNL);
        }
        integral += eL->weight(qL) * (_D[k] * inner_int_x * _base::template dNd< Proj>(eL, iL, qL, shiftL) +
                                      _D[2] * inner_int_y * _base::template dNd<!Proj>(eL, iL, qL, shiftL));
    }
    return integral;
}

template<class T, class I>
void structural_solver<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                                     Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                                     const T p1, const std::vector<bool>& inner_nodes) const {
    std::vector<std::unordered_set<I>> inner_portrait(K.rows()),
                                       bound_portrait(K_bound.rows());

    const auto indexator = [&inner_nodes, &inner_portrait, &bound_portrait, shift = 2 * first_node()](const I row, const I col) {
        if (inner_nodes[row] && inner_nodes[col]) {
            if (row <= col)
                inner_portrait[row - shift].insert(col);
        } else if (row != col) {
            if (!inner_nodes[col])
                bound_portrait[row - shift].insert(col);
        } else
            inner_portrait[row - shift].insert(col);
    };

    const auto structural_indexator = [&indexator](const I row, const I col) {
        indexator(row + I(X), col + I(X));
        indexator(row + I(X), col + I(Y));
        indexator(row + I(Y), col + I(X));
        indexator(row + I(Y), col + I(Y));
    };

    if (p1 > _base::MAX_LOCAL_WEIGHT) {
        _base::template mesh_run_loc(
            [this, &structural_indexator] (const size_t el, const size_t i, const size_t j) {
                structural_indexator(2 * mesh().node_number(el, i), 2 * mesh().node_number(el, j));
            });
    } else {
        _base::template mesh_run_nonloc(
            [this, &structural_indexator](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
                structural_indexator(2 * mesh().node_number(elL, iL), 2 * mesh().node_number(elNL, jNL));
            });
    }

    _base::convert_portrait(K, inner_portrait);
    _base::convert_portrait(K_bound, bound_portrait);
}

template<class T, class I>
template<class Influence_Function>
void structural_solver<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K,
                                          Eigen::SparseMatrix<T, Eigen::RowMajor, I>& K_bound,
                                          const T p1, const Influence_Function& influence_fun,
                                          const std::vector<bool>& inner_nodes) const {
    const auto filler_loc =
        [this, &K, &K_bound, &inner_nodes, p1, shift = 2 * first_node()]
                (const size_t el, const size_t i, const size_t j, const component proj, const component approx, const auto& integrate_rule) {
            const I row = 2 * mesh().node_number(el, i) + I(proj),
                    col = 2 * mesh().node_number(el, j) + I(approx);
            if (inner_nodes[row] && inner_nodes[col]) {
                if (row <= col)
                    K.coeffRef(row - shift, col) += p1 * integrate_rule(mesh().element_2d(el), i, j, quad_shift(el));
            } else if (row != col) {
                if (!inner_nodes[col])
                    K_bound.coeffRef(row - shift, col) += p1 * integrate_rule(mesh().element_2d(el), i, j, quad_shift(el));
            } else
                K.coeffRef(row - shift, col) = 1;
        };

    _base::template mesh_run_loc(
        [this, &filler_loc](const size_t el, const size_t i, const size_t j) {
#define SIGNATURE const Finite_Element_2D_Ptr& e, const size_t i, const size_t j, size_t quad_shift
            filler_loc(el, i, j, X, X, [this](SIGNATURE) { return integrate_loc<X, X>(e, i, j, quad_shift); });
            filler_loc(el, i, j, X, Y, [this](SIGNATURE) { return integrate_loc<X, Y>(e, i, j, quad_shift); });
            filler_loc(el, i, j, Y, X, [this](SIGNATURE) { return integrate_loc<Y, X>(e, i, j, quad_shift); });
            filler_loc(el, i, j, Y, Y, [this](SIGNATURE) { return integrate_loc<Y, Y>(e, i, j, quad_shift); });
#undef SIGNATURE
        });

    if (p1 < _base::MAX_LOCAL_WEIGHT) {
        const auto filler_nonloc =
            [this, &K, &K_bound, &inner_nodes, &influence_fun, p2 = 1 - p1, shift = 2 * first_node()]
                    (const size_t elL, const size_t iL, const size_t elNL, const size_t jNL, const component proj, const component approx, const auto& integrate_rule) {
                const I row = 2 * mesh().node_number(elL,  iL ) + I(proj),
                        col = 2 * mesh().node_number(elNL, jNL) + I(approx);
                if (inner_nodes[row] && inner_nodes[col]) {
                    if (row <= col)
                        K.coeffRef(row - shift, col) += p2 * integrate_rule(mesh().element_2d(elL ), iL,  quad_shift(elL ),
                                                                            mesh().element_2d(elNL), jNL, quad_shift(elNL), influence_fun);
                } else if (row != col)
                    if (!inner_nodes[col])
                        K_bound.coeffRef(row - shift, col) += p2 * integrate_rule(mesh().element_2d(elL ), iL,  quad_shift(elL ),
                                                                                  mesh().element_2d(elNL), jNL, quad_shift(elNL), influence_fun);
            };

        _base::template mesh_run_nonloc(
            [this, &filler_nonloc](const size_t elL, const size_t iL, const size_t elNL, const size_t jNL) {
#define SIGNATURE const Finite_Element_2D_Ptr& eL, const size_t iL, size_t shiftL, const Finite_Element_2D_Ptr& eNL, const size_t jNL, size_t shiftNL, const Influence_Function& influence_function
                filler_nonloc(elL, iL, elNL, jNL, X, X, [this](SIGNATURE) { return integrate_nonloc<X, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                filler_nonloc(elL, iL, elNL, jNL, X, Y, [this](SIGNATURE) { return integrate_nonloc<X, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                filler_nonloc(elL, iL, elNL, jNL, Y, X, [this](SIGNATURE) { return integrate_nonloc<Y, X, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
                filler_nonloc(elL, iL, elNL, jNL, Y, Y, [this](SIGNATURE) { return integrate_nonloc<Y, Y, Influence_Function>(eL, iL, shiftL, eNL, jNL, shiftNL, influence_function); });
#undef SIGNATURE
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
    create_matrix_portrait(K, K_bound, p1, inner_nodes);
    std::cout << "create_matrix_portrait: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    calc_matrix(K, K_bound, p1, influence_fun, inner_nodes);
    std::cout << "calc coeffs: " << omp_get_wtime() - time << std::endl;
}

template<class T, class I>
template<class Influence_Function>
solution<T, I> structural_solver<T, I>::stationary(const std::vector<bound_cond<T>> &bounds_cond,
                                                   const right_part<T>& right_part,
                                                   const T p1, const Influence_Function& influence_fun) {
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
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;
    //_base::PETSc_solver(f, K);

    return solution<T, I>{_base::get_mesh_info(), _D, u, p1, influence_fun};

    //std::cout << "x = " << std::endl;
    //VecView(x, nullptr);

//    time = omp_get_wtime();
//    //Eigen::PardisoLDLT<Eigen::SparseMatrix<T, Eigen::ColMajor, I>, Eigen::Lower> solver{K};
//    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
//    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
//    std::cout << "Matrix solve: " << omp_get_wtime() - time << std::endl;
//
//    return std::move(u);
}

template<class T, class I>
template<class Influence_Function>
void structural_solver<T, I>::temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                    const T alpha, const std::vector<T>& temperature,
                                                    const T p1, const Influence_Function& influence_fun) {
    const T c = alpha * _params.E / (1 - 2 * _params.nu);
    const auto [Tx, Ty] = _base::get_mesh_info()->calc_gradient(temperature);
    const std::vector<T> qTx = _base::get_mesh_info()->approx_in_quad(Tx),
                         qTy = _base::get_mesh_info()->approx_in_quad(Ty);

#pragma omp parallel for default(none) shared(f)
    for(size_t node = 0; node < mesh().nodes_count(); ++node) {
        T int_x = 0, int_y = 0;
        for(const I e : _base::nodes_elements_map(node)) {
            const auto& el = mesh().element_2d(e);
            const size_t i = _base::global_to_local_numbering(e).find(node)->second;
            for(size_t q = 0, shift = quad_shift(e); q < el->qnodes_count(); ++q, ++shift) {
                const T jac = el->weight(q) * jacobian(shift) * el->qN(i, q);
                int_x += qTx[shift] * jac;
                int_y += qTy[shift] * jac;
            }
        }
        f[2 * node]   += p1 * c * int_x;
        f[2 * node+1] += p1 * c * int_y;
    }

    if(p1 < _base::MAX_LOCAL_WEIGHT) {
        const T p2 = 1 - p1;
        for(size_t node = 0; node < mesh().nodes_count(); ++node) {
            T int_x = 0, int_y = 0;
            for(const I eL : _base::nodes_elements_map(node)) {
                const auto& elL = mesh().element_2d(eL);
                const size_t i = _base::global_to_local_numbering(eL).find(node)->second;
                for(const I eNL : _base::neighbors(eL)) {
                    const auto& elNL = mesh().element_2d(eNL);
                    for(size_t qL = 0, shiftL = quad_shift(eL); qL < elL->qnodes_count(); ++qL) {
                        T inner_int_x = 0, inner_int_y = 0;
                        for(size_t qNL = 0, shiftNL = quad_shift(eNL); qNL < elNL->qnodes_count(); ++qNL, ++shiftNL) {
                            const T influence_weight = elNL->weight(qNL) * jacobian(shiftNL) * influence_fun(quad_coord(shiftL), quad_coord(shiftNL));
                            inner_int_x += influence_weight * qTx[shiftNL];
                            inner_int_y += influence_weight * qTy[shiftNL];
                        }
                        const T jac = elL->weight(qL) * elL->qN(i, qL) * jacobian(shiftL);
                        int_x += jac * inner_int_x;
                        int_y += jac * inner_int_y;
                    }
                }
            }
            f[2 * node]   += p2 * c * int_x;
            f[2 * node+1] += p2 * c * int_y;
        }
    }
}

template<class T, class I>
template<class Influence_Function>
solution<T, I> structural_solver<T, I>::stationary(const std::vector<bound_cond<T>> &bounds_cond,
                                                   const right_part<T>& right_part,
                                                   const T alpha, const std::vector<T>& temperature,
                                                   const T p1, const Influence_Function& influence_fun) {
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
    temperature_condition(f, alpha, temperature, p1, influence_fun);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, I>, Eigen::Upper> solver{K};
    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

    return solution<T, I>{_base::get_mesh_info(), _D, u, p1, influence_fun};
}

}

#endif
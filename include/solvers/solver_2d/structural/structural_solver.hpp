#ifndef STRUCTURAL_SOLVER_HPP
#define STRUCTURAL_SOLVER_HPP

#include <functional>
#include <algorithm>
#include <omp.h>
#include "solver_2d/finite_element_solver_base_2d.hpp"
#include "structural_solution.hpp"

namespace nonlocal::structural {

enum class boundary_t : uint8_t {
    DISPLACEMENT = uint8_t(boundary_type::FIRST_KIND),
    PRESSURE     = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
using bound_cond = boundary_condition<T, boundary_t, 2>;

template<class T>
using right_part = right_partition<T, 2>;

template<class T, class I, class Matrix_Index>
class structural_solver : public finite_element_solver_base<T, I, Matrix_Index> {
    using _base = finite_element_solver_base<T, I, Matrix_Index>;
    using typename _base::component;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::jacobian;
    using _base::first_node;
    using _base::last_node;

    static void add_to_pair(std::array<T, 4>& pairs, const std::array<T, 2>& wdNd, const std::array<T, 2>& dNd) noexcept;
    static std::array<T, 4> calc_block(const std::array<T, 3>& D, const std::array<T, 4>& pairs) noexcept;

    std::array<T, 4> integrate_loc(const std::array<T, 3>& D, const size_t e, const size_t i, const size_t j) const;

    template<class Influence_Function>
    std::array<T, 4> integrate_nonloc(const std::array<T, 3>& D,
                                      const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                      const Influence_Function& influence_function) const;

    template<class Influence_Function>
    void temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                               const equation_parameters<T>& parameters,
                               const Influence_Function& influence_fun);

public:
    explicit structural_solver(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh);
    ~structural_solver() override = default;

    template<class Influence_Function>
    solution<T, I> stationary(const equation_parameters<T>& params, const std::vector<bound_cond<T>> &bounds_cond,
                              const right_part<T>& right_part, const Influence_Function& influence_fun);
};

template<class T, class I, class Matrix_Index>
structural_solver<T, I, Matrix_Index>::structural_solver(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class Matrix_Index>
void structural_solver<T, I, Matrix_Index>::add_to_pair(std::array<T, 4>& pairs, const std::array<T, 2>& wdNd, const std::array<T, 2>& dNd) noexcept {
    pairs[0] += wdNd[X] * dNd[X];
    pairs[1] += wdNd[X] * dNd[Y];
    pairs[2] += wdNd[Y] * dNd[X];
    pairs[3] += wdNd[Y] * dNd[Y];
}

template<class T, class I, class Matrix_Index>
std::array<T, 4> structural_solver<T, I, Matrix_Index>::calc_block(const std::array<T, 3>& D, const std::array<T, 4>& pairs) noexcept {
    return {
        D[0] * pairs[0] + D[2] * pairs[3],
        D[1] * pairs[1] + D[2] * pairs[2],
        D[1] * pairs[2] + D[2] * pairs[1],
        D[0] * pairs[3] + D[2] * pairs[0],
    };
}

template<class T, class I, class Matrix_Index>
std::array<T, 4> structural_solver<T, I, Matrix_Index>::integrate_loc(const std::array<T, 3>& D, const size_t e, const size_t i, const size_t j) const {
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    std::array<T, 4> pairs = {};
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj) {
        const T weight = el->weight(q) / jacobian(*J);
        const std::array<T, 2> wdNdi = {weight * (*dNdi)[X], weight * (*dNdi)[Y]};
        add_to_pair(pairs, wdNdi, *dNdj);
    }
    return calc_block(D, pairs);
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
std::array<T, 4> structural_solver<T, I, Matrix_Index>::integrate_nonloc(const std::array<T, 3>& D,
                                                                         const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                                                         const Influence_Function& influence_function) const {
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_start = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_start    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, 4> pairs = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_start;
        auto qcoordNL = qcoordNL_start;
        std::array<T, 2> inner_int = {};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_int[X] += influence_weight * (*dNdNL)[X];
            inner_int[Y] += influence_weight * (*dNdNL)[Y];
        }
        const std::array<T, 2> wdNdL = {elL->weight(qL) * (*dNdL)[X], elL->weight(qL) * (*dNdL)[Y]};
        add_to_pair(pairs, wdNdL, inner_int);
    }
    return calc_block(D, pairs);
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
void structural_solver<T, I, Matrix_Index>::temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                  const equation_parameters<T>& parameters,
                                                                  const Influence_Function& influence_fun) {
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
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
solution<T, I> structural_solver<T, I, Matrix_Index>::stationary(const equation_parameters<T>& parameters,
                                                                 const std::vector<bound_cond<T>> &bounds_cond,
                                                                 const right_part<T>& right_part,
                                                                 const Influence_Function& influence_fun) {
    double time = omp_get_wtime();
    const size_t rows = 2 * (last_node() - first_node()),
                 cols = 2 * mesh().nodes_count();
    const bool nonlocal_task = parameters.p1 < _base::MAX_LOCAL_WEIGHT;
    const std::vector<bool> inner_nodes = _base::template calc_inner_nodes(bounds_cond);
    const std::array<T, 3> D = hooke_matrix(parameters.nu, parameters.E, parameters.type);
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index> K_inner(rows, cols), K_bound(rows, cols);
    _base::template create_matrix_portrait<2>(K_inner, K_bound, inner_nodes, nonlocal_task);
    _base::template calc_matrix<2>(K_inner, K_bound, inner_nodes, nonlocal_task, influence_fun,
        [this, p1 = parameters.p1, &D](const size_t e, const size_t i, const size_t j) {
            using namespace metamath::function;
            return p1 * integrate_loc(D, e, i, j);
        },
        [this, p2 = 1 - parameters.p1, &D](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            using namespace metamath::function;
            return p2 * integrate_nonloc(D, eL, eNL, iL, jNL, influence_function);
        });

    std::cout << "Matrix create: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    if (parameters.thermoelasticity && parameters.delta_temperature.size() == mesh().nodes_count())
        temperature_condition(f, parameters, influence_fun);
    _base::template integrate_right_part(f, right_part);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> displacement = solver.solve(f);
    std::cout << "iterations = " << solver.iterations() << std::endl;
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

    return solution<T, I>{_base::mesh_proxy(), parameters, influence_fun, displacement};
}

}

#endif
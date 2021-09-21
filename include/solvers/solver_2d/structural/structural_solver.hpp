#ifndef STRUCTURAL_SOLVER_HPP
#define STRUCTURAL_SOLVER_HPP

#include "finite_element_solver_base_2d.hpp"
#include "structural_solution.hpp"
#include "conjugate_gradient.hpp"
#include <functional>
#include <algorithm>
#include <omp.h>

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
    using _base::X;
    using _base::Y;

    static void add_to_pair(std::array<T, 4>& pairs, const std::array<T, 2>& wdNd, const std::array<T, 2>& dNd) noexcept;
    static std::array<T, 4> calc_block(const std::array<T, 3>& D, const std::array<T, 4>& pairs) noexcept;

    std::array<T, 4> integrate_loc(const std::array<T, 3>& D, const size_t e, const size_t i, const size_t j) const;

    template<class Influence_Function>
    std::array<T, 4> integrate_nonloc(const std::array<T, 3>& D,
                                      const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                      const Influence_Function& influence_function) const;

    std::array<T, 2> integrate_temperature_loc(const std::vector<T>& eps_T, const size_t e, const size_t i) const;

    template<class Influence_Function>
    std::array<T, 2> integrate_temperature_nonloc(const std::vector<T>& eps_T, const size_t eL, const size_t eNL, const size_t iL,
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
        const T weight = el->weight(q) / _base::jacobian(*J);
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
    const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_begin    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, 4> pairs = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_begin;
        auto qcoordNL = qcoordNL_begin;
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
std::array<T, 2> structural_solver<T, I, Matrix_Index>::integrate_temperature_loc(const std::vector<T>& eps_T, const size_t e, const size_t i) const {
    std::array<T, 2> integral = {};
    const auto& el  = _base::mesh().element_2d(e);
          auto  dNd = _base::mesh_proxy()->dNdX(e, i);
          auto  eps = std::next(eps_T.cbegin(), _base::mesh_proxy()->quad_shift(e));
    for(size_t q = 0; q < el->qnodes_count(); ++q, ++dNd, ++eps) {
        const T weight = el->weight(q) * (*eps);
        integral[X] += weight * (*dNd)[X];
        integral[Y] += weight * (*dNd)[Y];
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
std::array<T, 2> structural_solver<T, I, Matrix_Index>::integrate_temperature_nonloc(const std::vector<T>& eps_T,
                                                                                     const size_t eL, const size_t eNL, const size_t iL,
                                                                                     const Influence_Function& influence_function) const {
    std::array<T, 2> integral = {};
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
    const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
    const auto  JNL_begin      = _base::mesh_proxy()->jacobi_matrix(eNL);
    const auto  epsNL_begin    = std::next(eps_T.cbegin(), _base::mesh_proxy()->quad_shift(eNL));
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto JNL = JNL_begin;
        auto epsNL = epsNL_begin;
        auto qcoordNL = qcoordNL_begin;
        T inner_integral = T{0};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++JNL, ++qcoordNL, ++epsNL)
            inner_integral += elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL) * (*epsNL) * _base::jacobian(*JNL);
        integral[X] += elL->weight(qL) * (*dNdL)[X] * inner_integral;
        integral[Y] += elL->weight(qL) * (*dNdL)[Y] * inner_integral;
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
void structural_solver<T, I, Matrix_Index>::temperature_condition(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                  const equation_parameters<T>& parameters,
                                                                  const Influence_Function& influence_fun) {
    const T nu = parameters.poisson(),
            E  = parameters.young();
    const T factor = 0.5 * parameters.alpha * E / (1 - nu);
    using namespace metamath::function;
    const std::vector<T> eps_T = factor * _base::mesh_proxy()->approx_in_quad(parameters.delta_temperature);

#pragma omp parallel for default(none) shared(f, eps_T, parameters)
    for(size_t node = _base::first_node(); node < _base::last_node(); ++node)
        for(const I e : _base::mesh_proxy()->nodes_elements_map(node)) {
            const size_t i = _base::mesh_proxy()->global_to_local_numbering(e).find(node)->second;
            const std::array<T, 2> integral = integrate_temperature_loc(eps_T, e, i);
            f[2 * node + X] += parameters.p1 * integral[X];
            f[2 * node + Y] += parameters.p1 * integral[Y];
        }

    if (parameters.p1 < _base::MAX_LOCAL_WEIGHT) {
        const size_t p2 = 1 - parameters.p1;
#pragma omp parallel for default(none) shared(f, eps_T, p2) firstprivate(influence_fun)
        for(size_t node = _base::first_node(); node < _base::last_node(); ++node)
            for(const I eL : _base::mesh_proxy()->nodes_elements_map(node)) {
                const size_t iL = _base::mesh_proxy()->global_to_local_numbering(eL).find(node)->second;
                for(const I eNL : _base::mesh_proxy()->neighbors(eL)) {
                    const std::array<T, 2> integral = integrate_temperature_nonloc(eps_T, eL, eNL, iL, influence_fun);
                    f[2 * node + X] += p2 * integral[X];
                    f[2 * node + Y] += p2 * integral[Y];
                }
            }
    }
}

template<class T, class I, class Matrix_Index>
template<class Influence_Function>
solution<T, I> structural_solver<T, I, Matrix_Index>::stationary(const equation_parameters<T>& parameters,
                                                                 const std::vector<bound_cond<T>> &bounds_cond,
                                                                 const right_part<T>& right_part,
                                                                 const Influence_Function& influence_fun) {
    double time = omp_get_wtime();
    const size_t rows = 2 * (_base::last_node() - _base::first_node()),
                 cols = 2 * _base::mesh().nodes_count();
    const bool nonlocal_task = parameters.p1 < _base::MAX_LOCAL_WEIGHT;
    const std::vector<bool> inner_nodes = _base::template calc_inner_nodes(bounds_cond);
    const std::array<T, 3> D = hooke_matrix(parameters);
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
    if (parameters.thermoelasticity && parameters.delta_temperature.size() == _base::mesh().nodes_count())
        temperature_condition(f, parameters, influence_fun);
    _base::template integrate_right_part(f, right_part);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
    std::cout << "Boundary cond: " << omp_get_wtime() - time << std::endl;

    time = omp_get_wtime();
    MPI_utils::MPI_ranges ranges = _base::mesh_proxy()->ranges();
    for(size_t rank = 0; rank < ranges.ranges().size(); ++rank)
        for(size_t& val : ranges.range(rank))
            val += val;
    const slae::conjugate_gradient<T, Matrix_Index> solver{K_inner, ranges};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> displacement = solver.solve(f);
    std::cout << "iterations = " << solver.iterations() << std::endl;
    std::cout << "residual = " << solver.residual() << std::endl;
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

    return solution<T, I>{_base::mesh_proxy(), parameters, influence_fun, displacement};
}

}

#endif
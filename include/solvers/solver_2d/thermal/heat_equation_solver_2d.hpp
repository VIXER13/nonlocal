#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <iostream>
#include <algorithm>
#include <omp.h>
#include "solver_2d/finite_element_solver_base_2d.hpp"
#include "heat_equation_parameters.hpp"
#include "heat_equation_solution_2d.hpp"

namespace nonlocal::heat {

enum class boundary_t : uint8_t {
    TEMPERATURE = uint8_t(boundary_type::FIRST_KIND),
    FLOW        = uint8_t(boundary_type::SECOND_KIND)
};

template<class T>
using bound_cond = boundary_condition<T, boundary_t, 1>;

template<class T>
using right_part = right_partition<T, 1>;

template<class T, class I, class Matrix_Index>
class heat_equation_solver_2d : public finite_element_solver_base<T, I, Matrix_Index> {
    using _base = finite_element_solver_base<T, I, Matrix_Index>;
    using _base::size;
    using _base::rank;
    using _base::first_node;
    using _base::last_node;
    using _base::X;
    using _base::Y;
    using _base::mesh;
    using _base::jacobian;

    std::array<T, 2> _lambda = { T{1}, T{1} };

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;
    template<calc_type Calc_Type>
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<calc_type Calc_Type, class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                const std::vector<bool>& inner_nodes,
                                const bool neumann_task, const bool nonlocal_task) const;

    void neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc,
                     const bool neumann_task,
                     const bool nonlocal_task, const Influence_Function& influence_fun,
                     const std::vector<bool>& inner_nodes) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                       Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const bool neumann_task,
                       const Integrate_Loc& integrate_rule_loc,
                       const Integrate_Nonloc& integrate_rule_nonloc,
                       const bool nonlocal_task, const Influence_Function& influence_fun) const;

    void nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                     const solver_parameters<T>& sol_parameters, const uintmax_t step);

public:
    explicit heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh);
    ~heat_equation_solver_2d() override = default;

    template<calc_type Calc_Type, class Influence_Function>
    solution<T, I> stationary(const equation_parameters<T, Calc_Type>& eq_parameters,
                              const std::vector<bound_cond<T>>& bounds_cond, const right_part<T>& right_part,
                              const T p1, const Influence_Function& influence_fun, const T volume = 0);

    template<calc_type Calc_Type, class Init_Distribution, class Right_Part, class Influence_Function>
    void nonstationary(const solver_parameters<T>& sol_parameters,
                       const equation_parameters<T, Calc_Type>& eq_parameters,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const Init_Distribution& init_dist, const Right_Part& right_part,
                       const T p1, const Influence_Function& influence_fun);
};

template<class T, class I, class Matrix_Index>
heat_equation_solver_2d<T, I, Matrix_Index>::heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh)
    : _base{mesh} {}

template<class T, class I, class Matrix_Index>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_basic(const size_t e, const size_t i) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q) * jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<calc_type Calc_Type>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    if constexpr (Calc_Type == calc_type::ISOTROPIC) {
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj)
            integral += el->weight(q) * ((*dNdi)[X] * (*dNdj)[X] + (*dNdi)[Y] * (*dNdj)[Y]) / jacobian(*J);
    } else if constexpr (Calc_Type == calc_type::ORTHOTROPIC) {
        std::array<T, 2> integral_part = {};
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj) {
            const T factor = el->weight(q) / jacobian(*J);
            integral_part[X] += factor * (*dNdi)[X] * (*dNdj)[X];
            integral_part[Y] += factor * (*dNdi)[Y] * (*dNdj)[Y];
        }
        integral = _lambda[X] * integral_part[X] + _lambda[Y] * integral_part[Y];
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<calc_type Calc_Type, class Influence_Function>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_nonloc(const size_t eL, const size_t eNL,
                                                             const size_t iL, const size_t jNL,
                                                             const Influence_Function& influence_function) const {
    T integral = 0;
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_start = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_start    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, Calc_Type == calc_type::ORTHOTROPIC ? 2 : 1> integral_part = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_start;
        auto qcoordNL = qcoordNL_start;
        std::array<T, 2> inner_integral_part = {};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_integral_part[X] += influence_weight * (*dNdNL)[X];
            inner_integral_part[Y] += influence_weight * (*dNdNL)[Y];
        }
        if constexpr (Calc_Type == calc_type::ISOTROPIC)
            integral += elL->weight(qL) * (inner_integral_part[X] * (*dNdL)[X] + inner_integral_part[Y] * (*dNdL)[Y]);
        else if constexpr (Calc_Type == calc_type::ORTHOTROPIC) {
            integral_part[X] += elL->weight(qL) * inner_integral_part[X] * (*dNdL)[X];
            integral_part[Y] += elL->weight(qL) * inner_integral_part[Y] * (*dNdL)[Y];
        }
    }
    if constexpr (Calc_Type == calc_type::ORTHOTROPIC)
        integral = _lambda[X] * integral_part[X] + _lambda[Y] * integral_part[Y];
    return integral;
}

template<class T, class I, class Matrix_Index>
void heat_equation_solver_2d<T, I, Matrix_Index>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                                         Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                                         const std::vector<bool>& inner_nodes,
                                                                         const bool neumann_task, const bool nonlocal_task) const {
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.outerIndexPtr()[row+1] = 1;
    _base::template create_matrix_portrait<1>(K_inner, K_bound, inner_nodes, nonlocal_task);
}

template<class T, class I, class Matrix_Index>
void heat_equation_solver_2d<T, I, Matrix_Index>::neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner) const {
#pragma omp parallel for default(none) shared(K_inner)
    for(size_t node = first_node(); node < last_node(); ++node) {
        T& val = K_inner.coeffRef(node - first_node(), mesh().nodes_count());
        for(const I e : _base::mesh_proxy()->nodes_elements_map(node))
            val += integrate_basic(e, _base::mesh_proxy()->global_to_local_numbering(e).find(node)->second);
    }
}

template<class T, class I, class Matrix_Index>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void heat_equation_solver_2d<T, I, Matrix_Index>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                           Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                           const Integrate_Loc& integrate_rule_loc,
                                                           const Integrate_Nonloc& integrate_rule_nonloc,
                                                           const bool neumann_task,
                                                           const bool nonlocal_task, const Influence_Function& influence_fun,
                                                           const std::vector<bool>& inner_nodes) const {
    _base::template mesh_run<_base::theory::LOCAL>(
        [this, &K_inner, &K_bound, &inner_nodes, &integrate_rule_loc](const size_t e, const size_t i, const size_t j) {
            const I row = mesh().node_number(e, i),
                    col = mesh().node_number(e, j);
            if (inner_nodes[row] && inner_nodes[col]) {
                if (row <= col)
                    K_inner.coeffRef(row - first_node(), col) += integrate_rule_loc(e, i, j);
            } else if (row != col) {
                if (!inner_nodes[col])
                    K_bound.coeffRef(row - first_node(), col) += integrate_rule_loc(e, i, j);
            } else
                K_inner.coeffRef(row - first_node(), col) = 1;
        }
    );

    if (nonlocal_task) {
        _base::template mesh_run<_base::theory::NONLOCAL>(
            [this, &K_inner, &K_bound, &inner_nodes, &integrate_rule_nonloc, &influence_fun]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const I row = mesh().node_number(eL,  iL ),
                        col = mesh().node_number(eNL, jNL);
                if (inner_nodes[row] && inner_nodes[col]) {
                    if (row <= col)
                        K_inner.coeffRef(row - first_node(), col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
                } else if (row != col)
                    if (!inner_nodes[col])
                        K_bound.coeffRef(row - first_node(), col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
            }
        );
    }
}

template<class T, class I, class Matrix_Index>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void heat_equation_solver_2d<T, I, Matrix_Index>::create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                             Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                             const std::vector<bound_cond<T>>& bounds_cond,
                                                             const bool neumann_task,
                                                             const Integrate_Loc& integrate_rule_loc,
                                                             const Integrate_Nonloc& integrate_rule_nonloc,
                                                             const bool nonlocal_task, const Influence_Function& influence_fun) const {
    std::vector<bool> inner_nodes(mesh().nodes_count(), true);
    _base::template boundary_nodes_run([this, &bounds_cond, &inner_nodes](const size_t b, const size_t el, const size_t i) {
        if(bounds_cond[b].type(0) == boundary_t::TEMPERATURE)
            inner_nodes[mesh().node_number(b, el, i)] = false;
    });

    double time = omp_get_wtime();
    create_matrix_portrait(K_inner, K_bound, inner_nodes, neumann_task, nonlocal_task);
    std::cout << "rank = " << rank() << std::endl;
    std::cout << "K.nonzero() = " << K_inner.nonZeros() << std::endl;
    std::cout << "K_bound.nonzero() = " << K_bound.nonZeros() << std::endl;
    std::cout << "create_matrix_portrait: " << omp_get_wtime() - time << std::endl << std::endl;

    time = omp_get_wtime();
    calc_matrix(K_inner, K_bound, integrate_rule_loc, integrate_rule_nonloc, neumann_task, nonlocal_task, influence_fun, inner_nodes);
    std::cout << "rank = " << rank() << std::endl;
    std::cout << "calc coeffs: " << omp_get_wtime() - time << std::endl;
}

template<class T, class I, class Matrix_Index>
template<calc_type Calc_Type, class Influence_Function>
solution<T, I> heat_equation_solver_2d<T, I, Matrix_Index>::stationary(const equation_parameters<T, Calc_Type>& eq_parameters,
                                                                    const std::vector<bound_cond<T>>& bounds_cond,
                                                                    const right_part<T>& right_part,
                                                                    const T p1, const Influence_Function& influence_fun,
                                                                    const T volume) {
    static constexpr auto boundary_checker = [](const bound_cond<T>& bound) { return bound.type(0) == boundary_t::FLOW; };
    const bool neumann_task = std::all_of(bounds_cond.cbegin(), bounds_cond.cend(), boundary_checker);
    const size_t rows = last_node() - first_node() + (neumann_task && rank() == size() - 1),
                 cols = mesh().nodes_count() + neumann_task;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    if(neumann_task) {
        //if(std::abs(std::accumulate(f.data(), f.data() + f.size(), T{0}, [](const T sum, const T val) { return sum + val; })) > 1e-5)
        //    throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
        f[mesh().nodes_count()] = volume;
    }

    if constexpr (Calc_Type == calc_type::ORTHOTROPIC)
        _lambda = eq_parameters.lambda;

    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index> K_inner(rows, cols), K_bound(rows, cols);
    const auto integrate_rule_loc = [this, factor_loc = Calc_Type == calc_type::ISOTROPIC ? p1 * eq_parameters.lambda[0] : p1]
                                    (const size_t e, const size_t i, const size_t j) {
        return factor_loc * integrate_loc<Calc_Type>(e, i, j);
    };
    const auto integrate_rule_nonloc =
        [this, factor_nonloc = Calc_Type == calc_type::ISOTROPIC ? (T{1} - p1) * eq_parameters.lambda[0] : (T{1} - p1)]
        (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            return factor_nonloc * integrate_nonloc<Calc_Type>(eL, eNL, iL, jNL, influence_function);
    };
    create_matrix(
        K_inner, K_bound, bounds_cond, neumann_task,
        integrate_rule_loc, integrate_rule_nonloc,
        p1 < _base::MAX_LOCAL_WEIGHT, influence_fun
    );
    if (neumann_task)
        neumann_task_col_fill(K_inner);

    _base::template integrate_right_part(f, right_part);
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);

    double time = omp_get_wtime();
//#ifdef MPI_USE
    //const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = _base::template MKL_solver<1>(f, K_inner);
    //Eigen::PardisoLLT<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{K_inner};
//    Eigen::Matrix<T, Eigen::Dynamic, 1> u = solver.solve(f);
//    f = u;
    //_base::PETSc_solver(f, K_inner);
//#else
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{K_inner};

    Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
//#endif
    std::cout << "System solve: " << omp_get_wtime() - time << std::endl;

    return solution<T, I>{_base::mesh_proxy(), temperature};
}

template<class T, class I, class Matrix_Index>
template<calc_type Calc_Type, class Init_Distribution, class Right_Part, class Influence_Function>
void heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary(const solver_parameters<T>& sol_parameters,
                                                             const equation_parameters<T, Calc_Type>& eq_parameters,
                                                             const std::vector<bound_cond<T>>& bounds_cond,
                                                             const Init_Distribution& init_dist, const Right_Part& right_part,
                                                             const T p1, const Influence_Function& influence_fun) {
    if constexpr (Calc_Type == calc_type::ORTHOTROPIC)
        _lambda = eq_parameters.lambda;

    static constexpr bool NOT_NEUMANN_TASK = false;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K_inner(mesh().nodes_count(), mesh().nodes_count()),
                                               K_bound(mesh().nodes_count(), mesh().nodes_count());
    const auto integrate_rule_loc = [this, factor_loc = Calc_Type == calc_type::ISOTROPIC ? p1 * eq_parameters.lambda[0] : p1]
                                    (const size_t e, const size_t i, const size_t j) {
        return factor_loc * integrate_loc<Calc_Type>(e, i, j);
    };
    const auto integrate_rule_nonloc =
        [this, factor_nonloc = Calc_Type == calc_type::ISOTROPIC ? (T{1} - p1) * eq_parameters.lambda[0] : (T{1} - p1)]
        (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            return factor_nonloc * integrate_nonloc<Calc_Type>(eL, eNL, iL, jNL, influence_function);
        };
    create_matrix(
        K_inner, K_bound, bounds_cond, NOT_NEUMANN_TASK,
        integrate_rule_loc, integrate_rule_nonloc,
        p1 < _base::MAX_LOCAL_WEIGHT, influence_fun
    );

    static constexpr bool LOCAL = false;
    Eigen::SparseMatrix<T, Eigen::RowMajor, I> C_inner(mesh().nodes_count(), mesh().nodes_count()),
                                               C_bound(mesh().nodes_count(), mesh().nodes_count());
    create_matrix(
        C_inner, C_bound, bounds_cond, NOT_NEUMANN_TASK,
        [this](const size_t e, const size_t i, const size_t j) { return integrate_basic_pair(e, i, j); },
        [](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) { return 0; },
        LOCAL, influence_fun
    );

    const T tau = (sol_parameters.time_interval.back() - sol_parameters.time_interval.front()) / sol_parameters.steps;
    C_bound.setZero();
    C_inner *= eq_parameters.rho * eq_parameters.c;
    K_bound *= tau;
    K_inner *= tau;
    K_inner += C_inner;
    _base::template boundary_nodes_run([this, &bounds_cond, &K_inner](const size_t b, const size_t el, const size_t i) {
        if(bounds_cond[b].type(0) == boundary_t::TEMPERATURE)
            K_inner.coeffRef(mesh().node_number(b, el, i), mesh().node_number(b, el, i)) = T{1};
    });

    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(mesh().nodes_count()),
                                        temperature_prev(mesh().nodes_count()),
                                        temperature_curr(mesh().nodes_count());
    for(size_t i = 0; i < mesh().nodes_count(); ++i)
        temperature_prev[i] = init_dist(mesh().node(i));

    Eigen::SparseMatrix<T, Eigen::ColMajor, Matrix_Index> KK = K_inner.transpose();
    K_inner.setZero();
    Eigen::SimplicialLLT<Eigen::SparseMatrix<T, Eigen::ColMajor, Matrix_Index>, Eigen::Lower> solver{KK};
    //Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{K_inner};
    if(sol_parameters.save_freq != std::numeric_limits<uintmax_t>::max())
        nonstationary_solver_logger(temperature_prev, sol_parameters, 0);
    for(size_t step = 1; step < sol_parameters.steps + 1; ++step) {
        f.setZero();
        _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
        _base::template integrate_right_part(f, right_partition<T, 1>{right_part});
        f *= tau;
        f += C_inner.template selfadjointView<Eigen::Upper>() * temperature_prev;
        _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);
        temperature_curr = solver.solve(f);
        temperature_prev.swap(temperature_curr);
        if(step % sol_parameters.save_freq == 0)
            nonstationary_solver_logger(temperature_prev, sol_parameters, step);
    }
}

template<class T, class I, class Matrix_Index>
void heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                                                           const solver_parameters<T>& sol_parameters, const uintmax_t step) {
    std::cout << "step = " << step << std::endl;
    if(sol_parameters.save_vtk)
        mesh::save_as_vtk(sol_parameters.save_path + std::to_string(step) + ".vtk", _base::mesh_proxy()->mesh(), temperature);
    if(sol_parameters.save_csv)
        mesh::save_as_csv(sol_parameters.save_path + std::to_string(step) + ".csv", _base::mesh_proxy()->mesh(), temperature);
    if(sol_parameters.calc_energy)
        std::cout << "Energy = " << _base::mesh_proxy()->integrate_solution(temperature) << std::endl;
}

}

#endif
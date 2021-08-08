#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include <iostream>
#include <algorithm>
#include <omp.h>
#include "solver_2d/finite_element_solver_base_2d.hpp"
#include "heat_equation_parameters.hpp"
#include "heat_equation_solution_2d.hpp"
#include "conjugate_gradient.hpp"

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
    using _base::X;
    using _base::Y;

    enum class matrix : bool {THERMAL_CONDUCTIVITY, HEAT_CAPACITY};

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;

    template<material_t Material>
    T integrate_loc([[maybe_unused]] const std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1>& lambda,
                    const size_t e, const size_t i, const size_t j) const;

    template<material_t Material, class Influence_Function>
    T integrate_nonloc([[maybe_unused]] const std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1>& lambda,
                       const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                const std::vector<bool>& inner_nodes,
                                const bool neumann_task, const bool nonlocal_task) const;

    void neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner) const;

    template<matrix Type, material_t Material, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                     Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                     const equation_parameters<T, Material>& eq_parameters,
                     const std::vector<bool>& inner_nodes,
                     const bool nonlocal_task, const Influence_Function& influence_fun) const;


    void prepare_nonstationary_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                      Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                      Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_inner,
                                      Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_bound,
                                      const std::vector<bool>& inner_nodes,
                                      const T rho, const T c, const T tau) const;

    template<class Init_Dist, class Right_Part>
    void nonstationary_calc(const solver_parameters<T>& sol_parameters,
                            const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                            const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                            const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_inner,
                            const std::vector<bound_cond<T>>& bounds_cond,
                            const Init_Dist& init_dist, const Right_Part& right_part) const;

    void nonstationary_solver_logger(const Eigen::Matrix<T, Eigen::Dynamic, 1>& temperature,
                                     const solver_parameters<T>& sol_parameters, const uintmax_t step) const;

public:
    explicit heat_equation_solver_2d(const std::shared_ptr<mesh::mesh_proxy<T, I>>& mesh);
    ~heat_equation_solver_2d() override = default;

    template<material_t Material, class Right_Part, class Influence_Function>
    solution<T, I> stationary(const equation_parameters<T, Material>& eq_parameters,
                              const std::vector<bound_cond<T>>& bounds_cond, const Right_Part& right_part,
                              const Influence_Function& influence_fun);

    template<material_t Material, class Init_Distribution, class Right_Part, class Influence_Function>
    void nonstationary(const solver_parameters<T>& sol_parameters,
                       const equation_parameters<T, Material>& eq_parameters,
                       const std::vector<bound_cond<T>>& bounds_cond,
                       const Init_Distribution& init_dist, const Right_Part& right_part,
                       const Influence_Function& influence_fun);
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
        integral += el->weight(q) * el->qN(i, q) * _base::jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el = _base::mesh().element_2d(e);
          auto  J  = _base::mesh_proxy()->jacobi_matrix(e);
    for(size_t q = 0; q < el->nodes_count(); ++q, ++J)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q) * _base::jacobian(*J);
    return integral;
}

template<class T, class I, class Matrix_Index>
template<material_t Material>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_loc([[maybe_unused]] const std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1>& lambda,
                                                             const size_t e, const size_t i, const size_t j) const {
    T integral = 0;
    const auto& el   = _base::mesh().element_2d(e);
          auto  J    = _base::mesh_proxy()->jacobi_matrix(e);
          auto  dNdi = _base::mesh_proxy()->dNdX(e, i),
                dNdj = _base::mesh_proxy()->dNdX(e, j);
    if constexpr (Material == material_t::ISOTROPIC) {
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj)
            integral += el->weight(q) * ((*dNdi)[X] * (*dNdj)[X] + (*dNdi)[Y] * (*dNdj)[Y]) / _base::jacobian(*J);
    } else if constexpr (Material == material_t::ORTHOTROPIC) {
        std::array<T, 2> integral_part = {};
        for(size_t q = 0; q < el->qnodes_count(); ++q, ++J, ++dNdi, ++dNdj) {
            const T factor = el->weight(q) / _base::jacobian(*J);
            integral_part[X] += factor * (*dNdi)[X] * (*dNdj)[X];
            integral_part[Y] += factor * (*dNdi)[Y] * (*dNdj)[Y];
        }
        integral = lambda[X] * integral_part[X] + lambda[Y] * integral_part[Y];
    }
    return integral;
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Influence_Function>
T heat_equation_solver_2d<T, I, Matrix_Index>::integrate_nonloc([[maybe_unused]] const std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1>& lambda,
                                                                const size_t eL, const size_t eNL, const size_t iL, const size_t jNL,
                                                                const Influence_Function& influence_function) const {
    T integral = 0;
    const auto& elL            = _base::mesh().element_2d(eL ),
              & elNL           = _base::mesh().element_2d(eNL);
          auto  qcoordL        = _base::mesh_proxy()->quad_coord(eL);
          auto  dNdL           = _base::mesh_proxy()->dNdX(eL, iL);
    const auto  qcoordNL_begin = _base::mesh_proxy()->quad_coord(eNL);
    const auto  dNdNL_begin    = _base::mesh_proxy()->dNdX(eNL, jNL);
    std::array<T, Material == material_t::ORTHOTROPIC ? 2 : 1> integral_part = {};
    for(size_t qL = 0; qL < elL->qnodes_count(); ++qL, ++qcoordL, ++dNdL) {
        auto dNdNL    = dNdNL_begin;
        auto qcoordNL = qcoordNL_begin;
        std::array<T, 2> inner_integral_part = {};
        for(size_t qNL = 0; qNL < elNL->qnodes_count(); ++qNL, ++qcoordNL, ++dNdNL) {
            const T influence_weight = elNL->weight(qNL) * influence_function(*qcoordL, *qcoordNL);
            inner_integral_part[X] += influence_weight * (*dNdNL)[X];
            inner_integral_part[Y] += influence_weight * (*dNdNL)[Y];
        }
        if constexpr (Material == material_t::ISOTROPIC)
            integral += elL->weight(qL) * (inner_integral_part[X] * (*dNdL)[X] + inner_integral_part[Y] * (*dNdL)[Y]);
        else if constexpr (Material == material_t::ORTHOTROPIC) {
            integral_part[X] += elL->weight(qL) * inner_integral_part[X] * (*dNdL)[X];
            integral_part[Y] += elL->weight(qL) * inner_integral_part[Y] * (*dNdL)[Y];
        }
    }
    if constexpr (Material == material_t::ORTHOTROPIC)
        integral = lambda[X] * integral_part[X] + lambda[Y] * integral_part[Y];
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
    for(size_t node = _base::first_node(); node < _base::last_node(); ++node) {
        T& val = K_inner.coeffRef(node - _base::first_node(), _base::mesh().nodes_count());
        for(const I e : _base::mesh_proxy()->nodes_elements_map(node))
            val += integrate_basic(e, _base::mesh_proxy()->global_to_local_numbering(e).find(node)->second);
    }
}

template<class T, class I, class Matrix_Index>
template<typename heat_equation_solver_2d<T, I, Matrix_Index>::matrix Type, material_t Material, class Influence_Function>
void heat_equation_solver_2d<T, I, Matrix_Index>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                              Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                              const equation_parameters<T, Material>& eq_parameters,
                                                              const std::vector<bool>& inner_nodes,
                                                              const bool nonlocal_task, const Influence_Function& influence_fun) const {
    if constexpr (Type == matrix::THERMAL_CONDUCTIVITY) {
        const T factor_loc    = Material == material_t::ISOTROPIC ? eq_parameters.p1 * eq_parameters.lambda[0] : eq_parameters.p1,
                factor_nonloc = Material == material_t::ISOTROPIC ? (T{1} - eq_parameters.p1) * eq_parameters.lambda[0] : (T{1} - eq_parameters.p1);
        _base::template calc_matrix<1>(K_inner, K_bound, inner_nodes, nonlocal_task, influence_fun,
            [this, factor_loc, &lambda = eq_parameters.lambda]
            (const size_t e, const size_t i, const size_t j) { return factor_loc * integrate_loc<Material>(lambda, e, i, j); },
            [this, factor_nonloc, &lambda = eq_parameters.lambda]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
                return factor_nonloc * integrate_nonloc<Material>(lambda, eL, eNL, iL, jNL, influence_function);
        });
    } else if constexpr (Type == matrix::HEAT_CAPACITY)
        _base::template calc_matrix<1>(K_inner, K_bound, inner_nodes, nonlocal_task, influence_fun,
            [this](const size_t e, const size_t i, const size_t j) { return integrate_basic_pair(e, i, j); },
            [this](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) { return T{0}; });
}

template<class T, class I, class Matrix_Index>
void heat_equation_solver_2d<T, I, Matrix_Index>::prepare_nonstationary_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                                               Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                                               Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_inner,
                                                                               Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_bound,
                                                                               const std::vector<bool>& inner_nodes,
                                                                               const T rho, const T c, const T tau) const {
    C_bound.setZero();
    C_inner *= rho * c;
    K_bound *= tau;
    K_inner *= tau;
    K_inner += C_inner;
    for(size_t i = 0; i < inner_nodes.size(); ++i)
        if (!inner_nodes[i])
            K_inner.coeffRef(i, i) = T{1};
}

template<class T, class I, class Matrix_Index>
template<class Init_Dist, class Right_Part>
void heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary_calc(const solver_parameters<T>& sol_parameters,
                                                                     const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_inner,
                                                                     const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& K_bound,
                                                                     const Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>& C_inner,
                                                                     const std::vector<bound_cond<T>>& bounds_cond,
                                                                     const Init_Dist& init_dist, const Right_Part& right_part) const {
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(_base::mesh().nodes_count()),
                                        temperature_prev(_base::mesh().nodes_count()),
                                        temperature_curr(_base::mesh().nodes_count());
    for(size_t i = 0; i < _base::mesh().nodes_count(); ++i)
        temperature_prev[i] = init_dist(_base::mesh().node(i));
    Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index>, Eigen::Upper> solver{K_inner};
    if(sol_parameters.save_freq != std::numeric_limits<uintmax_t>::max())
        nonstationary_solver_logger(temperature_prev, sol_parameters, 0);
    const T tau = (sol_parameters.time_interval.back() - sol_parameters.time_interval.front()) / sol_parameters.steps;
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
                                                                              const solver_parameters<T>& sol_parameters, const uintmax_t step) const {
    std::cout << "step = " << step << std::endl;
    if(sol_parameters.save_vtk)
        mesh::save_as_vtk(sol_parameters.save_path + std::to_string(step) + ".vtk", _base::mesh_proxy()->mesh(), temperature);
    if(sol_parameters.save_csv)
        mesh::save_as_csv(sol_parameters.save_path + std::to_string(step) + ".csv", _base::mesh_proxy()->mesh(), temperature);
    if(sol_parameters.calc_energy)
        std::cout << "Energy = " << _base::mesh_proxy()->integrate_solution(temperature) << std::endl;
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Right_Part, class Influence_Function>
solution<T, I> heat_equation_solver_2d<T, I, Matrix_Index>::stationary(const equation_parameters<T, Material>& eq_parameters,
                                                                    const std::vector<bound_cond<T>>& bounds_cond,
                                                                    const Right_Part& right_part,
                                                                    const Influence_Function& influence_fun) {
    static constexpr auto boundary_checker = [](const bound_cond<T>& bound) noexcept { return bound.type(0) == boundary_t::FLOW; };
    const bool neumann_task = std::all_of(bounds_cond.cbegin(), bounds_cond.cend(), boundary_checker);
    const size_t rows = _base::last_node() - _base::first_node() + (neumann_task && MPI_utils::MPI_rank() == MPI_utils::MPI_size() - 1),
                 cols = _base::mesh().nodes_count() + neumann_task;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(rows);
    _base::template integrate_boundary_condition_second_kind(f, bounds_cond);
    if(neumann_task) {
        //if(std::abs(std::accumulate(f.data(), f.data() + f.size(), T{0}, [](const T sum, const T val) { return sum + val; })) > 1e-5)
        //    throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
        f[_base::mesh().nodes_count()] = eq_parameters.integral;
    }

    const bool nonlocal_task = eq_parameters.p1 < _base::MAX_LOCAL_WEIGHT;
    const std::vector<bool> inner_nodes = _base::template calc_inner_nodes(bounds_cond);
    Eigen::SparseMatrix<T, Eigen::RowMajor, Matrix_Index> K_inner(rows, cols), K_bound(rows, cols);
    create_matrix_portrait(K_inner, K_bound, inner_nodes, neumann_task, nonlocal_task);
    calc_matrix<matrix::THERMAL_CONDUCTIVITY>(K_inner, K_bound, eq_parameters, inner_nodes, nonlocal_task, influence_fun);
    if (neumann_task)
        neumann_task_col_fill(K_inner);

    _base::template integrate_right_part(f, right_partition<T, 1>{right_part});
    _base::template boundary_condition_first_kind(f, bounds_cond, K_bound);

    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = slae::conjugate_gradient(K_inner, f, {}, _base::mesh_proxy()->ranges());
    return solution<T, I>{_base::mesh_proxy(), temperature};
}

template<class T, class I, class Matrix_Index>
template<material_t Material, class Init_Distribution, class Right_Part, class Influence_Function>
void heat_equation_solver_2d<T, I, Matrix_Index>::nonstationary(const solver_parameters<T>& sol_parameters,
                                                             const equation_parameters<T, Material>& eq_parameters,
                                                             const std::vector<bound_cond<T>>& bounds_cond,
                                                             const Init_Distribution& init_dist, const Right_Part& right_part,
                                                             const Influence_Function& influence_fun) {
    static constexpr bool LOCAL = false;
    static constexpr bool NOT_NEUMANN_TASK = false;
    const size_t size = _base::mesh().nodes_count();
    const bool nonlocal_task = eq_parameters.p1 < _base::MAX_LOCAL_WEIGHT;
    const std::vector<bool> inner_nodes = _base::template calc_inner_nodes(bounds_cond);

    Eigen::SparseMatrix<T, Eigen::RowMajor, I> K_inner(size, size), K_bound(size, size);
    create_matrix_portrait(K_inner, K_bound, inner_nodes, NOT_NEUMANN_TASK, nonlocal_task);
    calc_matrix<matrix::THERMAL_CONDUCTIVITY>(K_inner, K_bound, eq_parameters, inner_nodes, nonlocal_task, influence_fun);

    Eigen::SparseMatrix<T, Eigen::RowMajor, I> C_inner(size, size), C_bound(size, size);
    create_matrix_portrait(C_inner, C_bound, inner_nodes, NOT_NEUMANN_TASK, LOCAL);
    calc_matrix<matrix::HEAT_CAPACITY>(C_inner, C_bound, eq_parameters, inner_nodes, LOCAL, influence_fun);

    const T tau = (sol_parameters.time_interval.back() - sol_parameters.time_interval.front()) / sol_parameters.steps;
    prepare_nonstationary_matrix(K_inner, K_bound, C_inner, C_bound, inner_nodes, eq_parameters.rho, eq_parameters.c, tau);
    nonstationary_calc(sol_parameters, K_inner, K_bound, C_inner, bounds_cond, init_dist, right_part);
}

}

#endif
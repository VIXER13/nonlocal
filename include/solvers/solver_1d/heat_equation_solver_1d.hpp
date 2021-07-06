#ifndef HEAT_EQUATION_SOLVER_HPP
#define HEAT_EQUATION_SOLVER_HPP

#include "finite_element_solver_base_1d.hpp"

namespace nonlocal::heat {

template<class T, class I>
class heat_equation_solver_1d : public finite_element_solver_base_1d<T, I> {
    enum class calc_type : bool {STATIONARY, NONSTATIONARY};

    using _base = finite_element_solver_base_1d<T, I>;
    using typename _base::stationary_boundary;
    using typename _base::nonstatinary_boundary;
    using _base::mesh;

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_basic_pair(const size_t e, const size_t i, const size_t j) const;
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                const stationary_boundary& bound_cond,
                                const bool neumann_task, const bool nonlocal_task) const;

    void neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner) const;

    template<calc_type type, class Influence_Function>
    void calc_matrix(const equation_parameters<T>& parameters,
                     Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                     std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                     const stationary_boundary& bound_cond,
                     const bool nonlocal_task, const Influence_Function& influence_function) const;

public:
    explicit heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh);
    ~heat_equation_solver_1d() override = default;

    template<class Right_Part, class Influence_Function>
    std::vector<T> stationary(const equation_parameters<T>& parameters,
                              const stationary_boundary& bound_cond,
                              const Right_Part& right_part,
                              const Influence_Function& influence_function) const;
};

template<class T, class I>
heat_equation_solver_1d<T, I>::heat_equation_solver_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh)
    : finite_element_solver_base_1d<T, I>{mesh} {}

template<class T, class I>
T heat_equation_solver_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q);
    return integral * mesh()->jacobian();
}

template<class T, class I>
T heat_equation_solver_1d<T, I>::integrate_basic_pair(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->nodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * el->qN(j, q);
    return integral * mesh()->jacobian();
}

template<class T, class I>
T heat_equation_solver_1d<T, I>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qNxi(i, q) * el->qNxi(j, q);
    return integral / mesh()->jacobian();
}

template<class T, class I>
template<class Influence_Function>
T heat_equation_solver_1d<T, I>::integrate_nonloc(const size_t eL, const size_t eNL,
                                                  const size_t iL, const size_t jNL,
                                                  const Influence_Function& influence_function) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t qL = 0; qL < el->qnodes_count(); ++qL) {
        T inner_integral = T{0};
        const T qcoordL = mesh()->quad_coord(eL, qL);
        for(size_t qNL = 0; qNL < el->qnodes_count(); ++qNL) {
            const T qcoordNL = mesh()->quad_coord(eNL, qNL);
            inner_integral += el->weight(qNL) * influence_function(qcoordL, qcoordNL) * el->qNxi(jNL, qNL);
        }
        integral += el->weight(qL) * el->qNxi(iL, qL) * inner_integral;
    }
    return integral;
}

template<class T, class I>
void heat_equation_solver_1d<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                           const stationary_boundary& bound_cond,
                                                           const bool neumann_task, const bool nonlocal_task) const {
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.outerIndexPtr()[row+1] = 1;
    _base::create_matrix_portrait(K_inner, bound_cond, nonlocal_task);
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.innerIndexPtr()[K_inner.outerIndexPtr()[row+1]-1] = mesh()->nodes_count();
}

template<class T, class I>
void heat_equation_solver_1d<T, I>::neumann_task_col_fill(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner) const {
#pragma omp parallel for default(none) shared(K_inner)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node) {
        T& val = K_inner.coeffRef(node, mesh()->nodes_count());
        for(const auto& [e, i] : mesh()->node_elements(node).arr)
            if(e != std::numeric_limits<size_t>::max())
                val += integrate_basic(e, i);
    }
}

template<class T, class I>
template<typename heat_equation_solver_1d<T, I>::calc_type type, class Influence_Function>
void heat_equation_solver_1d<T, I>::calc_matrix(const equation_parameters<T>& parameters,
                                                Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                const stationary_boundary& bound_cond,
                                                const bool nonlocal_task, const Influence_Function& influence_function) const {
    if constexpr (type == calc_type::STATIONARY)
        _base::template calc_matrix(K_inner, K_bound, bound_cond, nonlocal_task, influence_function,
            [this, factor = parameters.lambda * parameters.p1](const size_t e, const size_t i, const size_t j) {
                return factor * integrate_loc(e, i, j);
            },
            [this, factor = parameters.lambda * (T{1} - parameters.p1)]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
                return factor * integrate_nonloc(eL, eNL, iL, jNL, influence_function);
            });
    else if constexpr (type == calc_type::NONSTATIONARY)
        _base::template create_matrix(K_inner, K_bound, bound_cond, nonlocal_task, influence_function,
            [this](const size_t e, const size_t i, const size_t j) { return integrate_basic_pair(e, i, j); },
            [](const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) { return 0; });
}

template<class T, class I>
template<class Right_Part, class Influence_Function>
std::vector<T> heat_equation_solver_1d<T, I>::stationary(const equation_parameters<T>& parameters,
                                                         const stationary_boundary& bound_cond,
                                                         const Right_Part& right_part,
                                                         const Influence_Function& influence_function) const {
    const bool neumann_task = bound_cond[0].first == boundary_condition_t::SECOND_KIND &&
                              bound_cond[1].first == boundary_condition_t::SECOND_KIND;
    if (neumann_task && bound_cond[0].second + bound_cond[1].second > 1e-5)
        throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};

    const size_t size = mesh()->nodes_count() + neumann_task;
    Eigen::SparseMatrix<T, Eigen::RowMajor> K_inner(size, size);
    std::array<std::vector<std::pair<size_t, T>>, 2> K_bound;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(size);

    const bool nonlocal_task = parameters.p1 < 0.999;
    create_matrix_portrait(K_inner, bound_cond, neumann_task, nonlocal_task);
    calc_matrix<calc_type::STATIONARY>(parameters, K_inner, K_bound, bound_cond, nonlocal_task, influence_function);
    if (neumann_task)
        neumann_task_col_fill(K_inner);

    _base::template integrate_right_part(f, right_part);
    _base::boundary_condition_second_kind(f, bound_cond);
    _base::boundary_condition_first_kind(f, bound_cond, K_bound);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh()->nodes_count())};
}

}

#endif
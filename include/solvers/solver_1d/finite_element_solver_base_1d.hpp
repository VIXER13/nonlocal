#ifndef FINITE_ELEMENT_SOLVER_BASE_1D_HPP
#define FINITE_ELEMENT_SOLVER_BASE_1D_HPP

#include "mesh.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

namespace nonlocal {

enum class boundary_condition_t : uint8_t {
    FIRST_KIND,
    SECOND_KIND
};

template<class T>
struct equation_parameters final {
    T lambda = T{1},
      rho    = T{1},
      c      = T{1},
      p1     = T{1},
      r      = T{0};
};

template<class T, class I>
class finite_element_solver_base_1d {
    std::shared_ptr<mesh::mesh_1d<T, I>> _mesh;

    enum class theory : bool { LOCAL, NONLOCAL };

    static void prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor>& K);

    T integrate_basic(const size_t e, const size_t i) const;
    T integrate_loc(const size_t e, const size_t i, const size_t j) const;
    template<class Influence_Function>
    T integrate_nonloc(const size_t eL, const size_t eNL,
                       const size_t iL, const size_t jNL,
                       const Influence_Function& influence_function) const;

    template<theory Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                const std::array<bool, 2> boundary_first_kind, const bool nonlocal_task) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                     std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                     const std::array<bool, 2> boundary_first_kind,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc,
                     const bool nonlocal_task, const Influence_Function& influence_fun) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                       std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                       const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                       const Integrate_Loc& integrate_rule_loc,
                       const Integrate_Nonloc& integrate_rule_nonloc,
                       const bool nonlocal_task, const Influence_Function& influence_fun) const;

    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                       const std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound) const;

    void boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond) const;

    template<class Function>
    T integrate_function(const size_t e, const size_t i, const Function& func) const;

    template<class Right_Part>
    void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const Right_Part& right_part) const;

public:
    explicit finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh);

    const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh() const;

    template<class Right_Part, class Influence_Function>
    std::vector<T> stationary(const equation_parameters<T>& parameters,
                              const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                              const Right_Part& right_part,
                              const Influence_Function& influence_function) const;
};

template<class T, class I>
void finite_element_solver_base_1d<T, I>::prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor>& K) {
    for(size_t i = 0; i < K.rows(); ++i)
        K.outerIndexPtr()[i+1] += K.outerIndexPtr()[i];
    K.data().resize(K.outerIndexPtr()[K.rows()]);
}

template<class T, class I>
finite_element_solver_base_1d<T, I>::finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T, I>>& finite_element_solver_base_1d<T, I>::mesh() const { return _mesh; }

template<class T, class I>
T finite_element_solver_base_1d<T, I>::integrate_basic(const size_t e, const size_t i) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q);
    return integral * mesh()->jacobian();
}

template<class T, class I>
T finite_element_solver_base_1d<T, I>::integrate_loc(const size_t e, const size_t i, const size_t j) const {
    T integral = T{0};
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qNxi(i, q) * el->qNxi(j, q);
    return integral / mesh()->jacobian();
}

template<class T, class I>
template<class Influence_Function>
T finite_element_solver_base_1d<T, I>::integrate_nonloc(const size_t eL, const size_t eNL,
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
template<typename finite_element_solver_base_1d<T, I>::theory Theory, class Callback>
void finite_element_solver_base_1d<T, I>::mesh_run(const Callback& callback) const {
#pragma omp parallel for default(none) firstprivate(callback) schedule(dynamic)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        for(const auto& [eL, iL] : mesh()->node_elements(node).arr)
            if(eL != std::numeric_limits<size_t>::max()) {
                if constexpr (Theory == theory::LOCAL)
                    for(size_t jL = 0; jL < mesh()->element()->nodes_count(); ++jL)
                        callback(eL, iL, jL);
                if constexpr (Theory == theory::NONLOCAL) {
                    const size_t finish = mesh()->right_neighbour(eL);
                    for(size_t eNL = mesh()->left_neighbour(eL); eNL < finish; ++eNL)
                        for(size_t jNL = 0; jNL < mesh()->element()->nodes_count(); ++jNL)
                            callback(eL, eNL, iL, jNL);
                }
            }
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                                 const std::array<bool, 2> boundary_first_kind, const bool nonlocal_task) const {
    for(size_t e = 0; e < mesh()->elements_count(); ++e) {
        const size_t right_neighbour = nonlocal_task ? mesh()->right_neighbour(e) : e + 1;
        for(size_t i = 0; i < mesh()->element()->nodes_count(); ++i) {
            const size_t row = e * (mesh()->element()->nodes_count() - 1) + i + 1;
            if (boundary_first_kind.front() && row == 1 || boundary_first_kind.back() && row == mesh()->nodes_count())
                K_inner.outerIndexPtr()[row] = 1;
            else {
                const bool last_node_first_kind = boundary_first_kind.back() && right_neighbour * (mesh()->element()->nodes_count() - 1) == mesh()->nodes_count()-1;
                K_inner.outerIndexPtr()[row] = (right_neighbour - e) * (mesh()->element()->nodes_count() - 1) - i + 1 - last_node_first_kind;
            }
        }
    }

    const bool neumann_task = !boundary_first_kind.front() && !boundary_first_kind.back();
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            ++K_inner.outerIndexPtr()[row+1];

    prepare_memory(K_inner);

    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            K_inner.innerIndexPtr()[K_inner.outerIndexPtr()[row+1]-1] = mesh()->nodes_count();

    for(size_t i = 0; i < K_inner.rows(); ++i)
        for(size_t j = K_inner.outerIndexPtr()[i], k = i; j < K_inner.outerIndexPtr()[i+1]; ++j, ++k) {
            if (j < K_inner.outerIndexPtr()[i+1] - neumann_task)
                K_inner.innerIndexPtr()[j] = k;
            K_inner.valuePtr()[j] = T{0};
        }
    //std::cout << Eigen::MatrixXd{K_inner} << std::endl << std::endl;
}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void finite_element_solver_base_1d<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                      std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                      const std::array<bool, 2> boundary_first_kind,
                                                      const Integrate_Loc& integrate_rule_loc,
                                                      const Integrate_Nonloc& integrate_rule_nonloc,
                                                      const bool nonlocal_task, const Influence_Function& influence_fun) const {
    if (!boundary_first_kind.front() && !boundary_first_kind.back()) {
#pragma omp parallel for default(none) shared(K_inner)
        for(size_t node = 0; node < mesh()->nodes_count(); ++node) {
            T& val = K_inner.coeffRef(node, mesh()->nodes_count());
            for(const auto& [e, i] : mesh()->node_elements(node).arr)
                if(e != std::numeric_limits<size_t>::max())
                    val += integrate_basic(e, i);
        }
    }

    mesh_run<theory::LOCAL>(
        [this, &K_inner, &K_bound, boundary_first_kind, &integrate_rule_loc](const size_t e, const size_t i, const size_t j) {
            const I row = mesh()->node_number(e, i),
                    col = mesh()->node_number(e, j);
            if (row <= col)
                K_inner.coeffRef(row, col) += integrate_rule_loc(e, i, j);
//            if (inner_nodes[row] && inner_nodes[col]) {
//                if (row <= col)
//                    K_inner.coeffRef(row - first_node(), col) += integrate_rule_loc(e, i, j);
//            } else if (row != col) {
//                if (!inner_nodes[col])
//                    K_bound.coeffRef(row - first_node(), col) += integrate_rule_loc(e, i, j);
//            } else
//                K_inner.coeffRef(row - first_node(), col) = 1;
        }
    );

    if (nonlocal_task) {
        std::cout << "test" << std::endl;
        mesh_run<theory::NONLOCAL>(
            [this, &K_inner, &K_bound, boundary_first_kind, &integrate_rule_nonloc, &influence_fun]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const I row = mesh()->node_number(eL,  iL ),
                        col = mesh()->node_number(eNL, jNL);
                if (row <= col)
                    K_inner.coeffRef(row, col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
//                if (inner_nodes[row] && inner_nodes[col]) {
//                    if (row <= col)
//                        K_inner.coeffRef(row - first_node(), col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
//                } else if (row != col)
//                    if (!inner_nodes[col])
//                        K_bound.coeffRef(row - first_node(), col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
            }
        );
    }
}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void finite_element_solver_base_1d<T, I>::create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                        std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                        const Integrate_Loc& integrate_rule_loc,
                                                        const Integrate_Nonloc& integrate_rule_nonloc,
                                                        const bool nonlocal_task, const Influence_Function& influence_fun) const {
    const std::array<bool, 2> boundary_first_kind = {bound_cond.front().first == boundary_condition_t::FIRST_KIND,
                                                     bound_cond.back().first  == boundary_condition_t::FIRST_KIND};
    create_matrix_portrait(K_inner, boundary_first_kind, nonlocal_task);
    calc_matrix(K_inner, K_bound, boundary_first_kind, integrate_rule_loc, integrate_rule_nonloc, nonlocal_task, influence_fun);
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                                        const std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound) const {
    std::array<T*, 2> fval = {&f[0], &f[mesh()->nodes_count()-1]};
    for(size_t b = 0; b < 2; ++b)
        if(bound_cond[b].first == boundary_condition_t::FIRST_KIND) {
            for(const auto& [i, val] : K_bound[b])
                f[i] += val * bound_cond[b].second;
            *fval[b] = bound_cond[b].second;
        }
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                         const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond) const {
    std::array<T*, 2> fval = {&f[0], &f[mesh()->nodes_count()-1]};
    for(size_t b = 0; b < 2; ++b)
        if(bound_cond[b].first == boundary_condition_t::SECOND_KIND)
            *fval[b] += bound_cond[b].second;
}

template<class T, class I>
template<class Function>
T finite_element_solver_base_1d<T, I>::integrate_function(const size_t e, const size_t i, const Function& func) const {
    T integral = 0;
    const auto& el = mesh()->element();
    for(size_t q = 0; q < el->qnodes_count(); ++q)
        integral += el->weight(q) * el->qN(i, q) * func(mesh()->quad_coord(e, q));
    return integral * mesh()->jacobian();
}

template<class T, class I>
template<class Right_Part>
void finite_element_solver_base_1d<T, I>::integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const Right_Part& right_part) const {
#pragma omp parallel for default(none) shared(f, right_part)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        for(const auto& [e, i] : mesh()->node_elements(node).arr)
            if (e != std::numeric_limits<size_t>::max())
                f[node] += integrate_function(e, i, right_part);
}

template<class T, class I>
template<class Right_Part, class Influence_Function>
std::vector<T> finite_element_solver_base_1d<T, I>::stationary(const equation_parameters<T>& parameters,
                                                               const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                               const Right_Part& right_part,
                                                               const Influence_Function& influence_function) const {
    const bool neumann_task = bound_cond[0].first == boundary_condition_t::SECOND_KIND &&
                              bound_cond[1].first == boundary_condition_t::SECOND_KIND;
    if (neumann_task && bound_cond[0].second + bound_cond[1].second > 1e-5)
        throw std::domain_error{"The problem is unsolvable. Contour integral != 0."};
    const bool nonlocal_task = parameters.p1 < 0.999;

    const size_t size = mesh()->nodes_count() + neumann_task;
    Eigen::SparseMatrix<T, Eigen::RowMajor> K_inner(size, size);
    std::array<std::vector<std::pair<size_t, T>>, 2> K_bound;
    Eigen::Matrix<T, Eigen::Dynamic, 1> f = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(size);

    const auto integrate_rule_loc = [this, factor = parameters.lambda * parameters.p1](const size_t e, const size_t i, const size_t j) {
        return factor * integrate_loc(e, i, j);
    };
    const auto integrate_rule_nonloc =
        [this, factor = parameters.lambda * (T{1} - parameters.p1)]
        (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL, const Influence_Function& influence_function) {
            return factor * integrate_nonloc(eL, eNL, iL, jNL, influence_function);
    };

    create_matrix(K_inner, K_bound, bound_cond, integrate_rule_loc, integrate_rule_nonloc, nonlocal_task, influence_function);
    //std::cout << Eigen::MatrixXd{K_inner} << std::endl << std::endl;

    integrate_right_part(f, right_part);
    boundary_condition_second_kind(f, bound_cond);
    boundary_condition_first_kind(f, bound_cond, K_bound);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh()->nodes_count())};
}

}

#endif
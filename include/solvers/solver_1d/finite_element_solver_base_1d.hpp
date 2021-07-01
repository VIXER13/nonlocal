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

    static void sort_indices(Eigen::SparseMatrix<T, Eigen::RowMajor>& K);
    static void prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor>& K);

    template<theory Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                const std::array<bool, 2> boundary_first_kind, const bool nonlocal_task) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                     std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc,
                     const bool neumann_task,
                     const bool nonlocal_task, const Influence_Function& influence_fun) const;

    template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
    void create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                       std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                       const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                       const bool neumann_task,
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
                                                                 std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                                 const std::array<bool, 2> boundary_first_kind, const bool nonlocal_task) const {
    for(size_t row = 0; row < K_inner.rows(); ++row)
        K_inner.outerIndexPtr()[row+1] = 1;

    const bool neumann_task = !boundary_first_kind.front() && !boundary_first_kind.back();
    if (neumann_task)
        for(size_t row = 0; row < K_inner.rows(); ++row)
            ++K_inner.outerIndexPtr()[row+1];

    std::array<size_t, 2> boundary_elements_count = {};
    mesh_run<theory::LOCAL>(
        [this, &K_inner, &boundary_elements_count, boundary_first_kind](const size_t e, const size_t i, const size_t j) {
            const size_t row = mesh()->node_number(e, i),
                         col = mesh()->node_number(e, j);
            if (boundary_first_kind.front() && (col == 0 || row == 0)) { // Left boundary
                if (row == 0)
                    ++boundary_elements_count[0];
            } else if (boundary_first_kind.back() && (col == mesh()->nodes_count()-1 || row == mesh()->nodes_count()-1)) { // Right boundary
                if (row == mesh()->nodes_count()-1)
                    ++boundary_elements_count[1];
            } else if (row < col) // Regular
                ++K_inner.outerIndexPtr()[row + 1];
        }
    );

    K_bound[0].reserve(boundary_elements_count[0]);
    K_bound[1].reserve(boundary_elements_count[1]);
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
}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void finite_element_solver_base_1d<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                      std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                      const Integrate_Loc& integrate_rule_loc,
                                                      const Integrate_Nonloc& integrate_rule_nonloc,
                                                      const bool neumann_task,
                                                      const bool nonlocal_task, const Influence_Function& influence_fun) const {

}

template<class T, class I>
template<class Integrate_Loc, class Integrate_Nonloc, class Influence_Function>
void finite_element_solver_base_1d<T, I>::create_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                        std::array<std::vector<std::pair<size_t, T>>, 2>& K_bound,
                                                        const std::array<std::pair<boundary_condition_t, T>, 2>& bound_cond,
                                                        const bool neumann_task,
                                                        const Integrate_Loc& integrate_rule_loc,
                                                        const Integrate_Nonloc& integrate_rule_nonloc,
                                                        const bool nonlocal_task, const Influence_Function& influence_fun) const {
    const std::array<bool, 2> boundary_first_kind = {bound_cond.front().first == boundary_condition_t::FIRST_KIND,
                                                     bound_cond.back().first  == boundary_condition_t::FIRST_KIND};
    create_matrix_portrait(K_inner, K_bound, boundary_first_kind, nonlocal_task);
    calc_matrix(K_inner, K_bound, integrate_rule_loc, integrate_rule_nonloc, neumann_task, nonlocal_task, influence_fun);
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

    create_matrix(K_inner, K_bound, bound_cond, neumann_task,
                  [](){}, [](){}, nonlocal_task, influence_function);

    boundary_condition_second_kind(f, bound_cond);
    integrate_right_part(f, right_part);
    boundary_condition_first_kind(f, bound_cond, K_bound);

    const Eigen::ConjugateGradient<Eigen::SparseMatrix<T, Eigen::RowMajor>, Eigen::Upper> solver{K_inner};
    const Eigen::Matrix<T, Eigen::Dynamic, 1> temperature = solver.solve(f);
    return std::vector<T>{temperature.cbegin(), std::next(temperature.cbegin(), mesh()->nodes_count())};
}

}

#endif
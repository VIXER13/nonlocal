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
struct solver_parameters final {
    std::string save_path; // Путь куда сохранять данные
    std::array<T, 2> time_interval = {0, 1};
    uintmax_t steps = 100,
              save_freq = 1; // Частота сохранения
    bool save_csv    = true, // Сохранять .csv файлы в формате (x1, x2, T)
         calc_energy = true; // Вычислять энергия при сохранении, иногда полезно для контроля расчёта
};

template<class T, class I>
class finite_element_solver_base_1d {
    std::shared_ptr<mesh::mesh_1d<T, I>> _mesh;

protected:
    enum class theory : bool { LOCAL, NONLOCAL };

    using stationary_boundary = std::array<std::pair<boundary_condition_t, T>, 2>;
    using nonstatinary_boundary = std::array<std::pair<boundary_condition_t, std::function<T(T)>>, 2>;

    explicit finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh);

    static void prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor>& K);
    static stationary_boundary convert_nonstationary_boundary_to_stationary(const nonstatinary_boundary& bound_cond, const T t);

    void create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                const stationary_boundary& bound_cond,
                                const bool nonlocal_task) const;

    template<theory Theory, class Callback>
    void mesh_run(const Callback& callback) const;

    template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
    void calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                     std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                     const stationary_boundary& bound_cond,
                     const bool nonlocal_task, const Influence_Function& influence_fun,
                     const Integrate_Loc& integrate_rule_loc,
                     const Integrate_Nonloc& integrate_rule_nonloc) const;

    void boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                       const stationary_boundary& bound_cond,
                                       const std::array<std::unordered_map<size_t, T>, 2>& K_bound) const;

    void boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                        const stationary_boundary& bound_cond) const;

    template<class Function>
    T integrate_function(const size_t e, const size_t i, const Function& func) const;

    template<class Right_Part>
    void integrate_right_part(Eigen::Matrix<T, Eigen::Dynamic, 1>& f, const Right_Part& right_part) const;

public:
    virtual ~finite_element_solver_base_1d() = default;

    const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh() const;
};

template<class T, class I>
void finite_element_solver_base_1d<T, I>::prepare_memory(Eigen::SparseMatrix<T, Eigen::RowMajor>& K) {
    for(size_t i = 0; i < K.rows(); ++i)
        K.outerIndexPtr()[i+1] += K.outerIndexPtr()[i];
    K.data().resize(K.outerIndexPtr()[K.rows()]);
    for(size_t i = 0; i < K.outerIndexPtr()[K.rows()]; ++i) {
        K.innerIndexPtr()[i] = std::numeric_limits<I>::max();
        K.valuePtr()[i] = T{0};
    }
}

template<class T, class I>
typename finite_element_solver_base_1d<T, I>::stationary_boundary
finite_element_solver_base_1d<T, I>::convert_nonstationary_boundary_to_stationary(const nonstatinary_boundary& bound_cond, const T t) {
    return {
        std::pair{bound_cond[0].first, bound_cond[0].second(t)},
        std::pair{bound_cond[1].first, bound_cond[1].second(t)},
    };
}

template<class T, class I>
finite_element_solver_base_1d<T, I>::finite_element_solver_base_1d(const std::shared_ptr<mesh::mesh_1d<T, I>>& mesh)
    : _mesh{mesh} {}

template<class T, class I>
const std::shared_ptr<mesh::mesh_1d<T, I>>& finite_element_solver_base_1d<T, I>::mesh() const { return _mesh; }

template<class T, class I>
void finite_element_solver_base_1d<T, I>::create_matrix_portrait(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                                 const stationary_boundary& bound_cond,
                                                                 const bool nonlocal_task) const {
#pragma omp parallel for default(none) shared(K_inner, bound_cond, nonlocal_task)
    for(size_t node = 0; node < mesh()->nodes_count(); ++node)
        if (bound_cond.front().first == boundary_condition_t::FIRST_KIND && node == 0 ||
            bound_cond.back ().first == boundary_condition_t::FIRST_KIND && node == mesh()->nodes_count()-1)
            K_inner.outerIndexPtr()[node+1] = 1;
        else {
            const auto [eL, iL, eR, iR] = mesh()->node_elements(node).named;
            const size_t e = eR == std::numeric_limits<size_t>::max() ? eL : eR,
                         i = iR == std::numeric_limits<size_t>::max() ? iL : iR,
                         right_neighbour = nonlocal_task ? mesh()->right_neighbour(e) : e + 1;
            const bool last_node_first_kind = bound_cond.back().first  == boundary_condition_t::FIRST_KIND &&
                                              right_neighbour * (mesh()->element()->nodes_count() - 1) == mesh()->nodes_count()-1;
            K_inner.outerIndexPtr()[node+1] += (right_neighbour - e) * (mesh()->element()->nodes_count() - 1) - i + 1 - last_node_first_kind;
        }

    prepare_memory(K_inner);

    for(size_t i = 0; i < K_inner.rows(); ++i)
        for(size_t j = K_inner.outerIndexPtr()[i], k = i; j < K_inner.outerIndexPtr()[i+1]; ++j, ++k)
            K_inner.innerIndexPtr()[j] = k;
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
template<class Influence_Function, class Integrate_Loc, class Integrate_Nonloc>
void finite_element_solver_base_1d<T, I>::calc_matrix(Eigen::SparseMatrix<T, Eigen::RowMajor>& K_inner,
                                                      std::array<std::unordered_map<size_t, T>, 2>& K_bound,
                                                      const stationary_boundary& bound_cond,
                                                      const bool nonlocal_task, const Influence_Function& influence_fun,
                                                      const Integrate_Loc& integrate_rule_loc,
                                                      const Integrate_Nonloc& integrate_rule_nonloc) const {
    const auto boundary_calc = [&K_bound, &K_inner](const size_t b, const size_t row, const size_t col, const T integral) {
        if (col == row)
            K_inner.coeffRef(row, col) = 1;
        else {
            const auto [it, flag] = K_bound[b].template try_emplace(col, integral);
            if (!flag)
                it->second += integral;
        }
    };

    const std::array<bool, 2> boundary_first_kind = {bound_cond.front().first == boundary_condition_t::FIRST_KIND,
                                                     bound_cond.back().first  == boundary_condition_t::FIRST_KIND};

    mesh_run<theory::LOCAL>(
        [this, &K_inner, &K_bound, boundary_first_kind, &integrate_rule_loc, &boundary_calc](const size_t e, const size_t i, const size_t j) {
            const I row = mesh()->node_number(e, i),
                    col = mesh()->node_number(e, j);
            if (boundary_first_kind.front() && (row == 0 || col == 0)) {
                if (row == 0)
                    boundary_calc(0, row, col, integrate_rule_loc(e, i, j));
            } else if (boundary_first_kind.back() && (row == mesh()->nodes_count()-1 || col == mesh()->nodes_count()-1)) {
                if (row == mesh()->nodes_count()-1)
                    boundary_calc(1, row, col, integrate_rule_loc(e, i, j));
            } else if (row <= col)
                K_inner.coeffRef(row, col) += integrate_rule_loc(e, i, j);
        }
    );

    if (nonlocal_task) {
        mesh_run<theory::NONLOCAL>(
            [this, &K_inner, &K_bound, boundary_first_kind, &boundary_calc, &integrate_rule_nonloc, &influence_fun]
            (const size_t eL, const size_t eNL, const size_t iL, const size_t jNL) {
                const I row = mesh()->node_number(eL,  iL ),
                        col = mesh()->node_number(eNL, jNL);
                if (boundary_first_kind.front() && (row == 0 || col == 0)) {
                    if (row == 0)
                        boundary_calc(0, row, col, integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun));
                } else if (boundary_first_kind.back() && (row == mesh()->nodes_count()-1 || col == mesh()->nodes_count()-1)) {
                    if (row == mesh()->nodes_count()-1)
                        boundary_calc(1, row, col, integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun));
                } else if (row <= col)
                    K_inner.coeffRef(row, col) += integrate_rule_nonloc(eL, eNL, iL, jNL, influence_fun);
            }
        );
    }
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_first_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                        const stationary_boundary& bound_cond,
                                                                        const std::array<std::unordered_map<size_t, T>, 2>& K_bound) const {
    std::array<T*, 2> fval = {&f[0], &f[mesh()->nodes_count()-1]};
    for(size_t b = 0; b < 2; ++b)
        if(bound_cond[b].first == boundary_condition_t::FIRST_KIND) {
            for(const auto& [i, val] : K_bound[b])
                f[i] -= val * bound_cond[b].second;
            *fval[b] = bound_cond[b].second;
        }
}

template<class T, class I>
void finite_element_solver_base_1d<T, I>::boundary_condition_second_kind(Eigen::Matrix<T, Eigen::Dynamic, 1>& f,
                                                                         const stationary_boundary& bound_cond) const {
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

}

#endif